function validate_dgp(A, B, cfg, varargin)
    % VALIDATE_DGP - quick consistency checks & assertions for the DGP outputs
    % Usage:
    %   validate_dgp(A,B,cfg);                 % warnings only
    %   validate_dgp(A,B,cfg,'strict',true);   % throw error on failure
    
    p = inputParser;
    addParameter(p,'strict',false,@islogical);
    parse(p,varargin{:});
    strict = p.Results.strict;
    
    Aok = A(~A.skipped,:);
    if isempty(Aok), warning('No non-skipped auctions.'); return; end
    
    % 1) Winner uniqueness per auction
    g = groupsummary(B(B.winner,:),'auction_id','sum','winner');
    bad = g.sum_winner ~= 1;
    if any(bad)
        msg = sprintf('[VALIDATE] Non-unique winner in %d auctions.', sum(bad));
        if strict, error(msg); else, warning(msg); end
    end
    
    % 2) Winner score equals min entrant score
    auctions = unique(Aok.auction_id);
    badCnt = 0;
    for k = 1:numel(auctions)
        aid = auctions(k);
        rows = B.auction_id==aid & B.entered;
        if ~any(rows), continue; end
        smin = min(B.score(rows));
        sA   = Aok.win_score(Aok.auction_id==aid);
        if abs(smin - sA) > 1e-8
            badCnt = badCnt + 1;
        end
    end
    if badCnt>0
        msg = sprintf('[VALIDATE] Winner score mismatch in %d auctions.', badCnt);
        if strict, error(msg); else, warning(msg); end
    end
    
    % 3) Non-entrants must have NaN bids
    mask = ~B.entered;
    if ~(all(isnan(B.score(mask))) && all(isnan(B.b0(mask))) && all(isnan(B.b2(mask))))
        msg = '[VALIDATE] Non-entrants should have NaN for bids.';
        if strict, error(msg); else, warning(msg); end
    end
    
    % 4) Threshold entry consistency: entered iff k <= kbar
    bad = false;
    for k = 1:height(Aok)
        aid = Aok.auction_id(k); thr = Aok.kbar(k);
        rows = B.auction_id==aid;
        bad = bad | any( (B.k(rows) <= thr) ~= B.entered(rows) );
    end
    if bad
        msg = '[VALIDATE] Entry threshold inconsistency (k vs kbar).';
        if strict, error(msg); else, warning(msg); end
    end
    
    % 5) Overrun & final payment identities for winners
    bad = false;
    for k = 1:height(Aok)
        aid = Aok.auction_id(k);
        rA  = Aok.auction_id==aid;
        rw  = B.auction_id==aid & B.winner;
        if ~any(rw), continue; end
        over_calc = B.b2(rw) .* ( (B.e2(rw) + Aok.epsilon(rA)) - 1 );
        pay_calc  = B.b0(rw) + B.b2(rw) .* (B.e2(rw) + Aok.epsilon(rA));
        bad = bad | abs(over_calc - Aok.overrun(rA)) > 1e-8 ...
                  | abs(pay_calc  - Aok.pay_final(rA)) > 1e-8;
    end
    if bad
        msg = '[VALIDATE] Winner identities (overrun/pay_final) violated.';
        if strict, error(msg); else, warning(msg); end
    end
    
    % 6) (Optional) monotonicity of s(c) per auction (non-decreasing in c)
    %     NOTE: This is a numerical property; allow small tolerance.
    badCnt = 0;
    for k = 1:height(Aok)
        aid  = Aok.auction_id(k);
        rows = B.auction_id==aid & B.entered;
        if ~any(rows), continue; end
        T = sortrows(B(rows,{'c','score'}),'c');
        if any(diff(T.score) < -1e-6)
            badCnt = badCnt + 1;
        end
    end
    if badCnt>0
        msg = sprintf('[VALIDATE] s(c) monotonicity failed in %d auctions (tolerance 1e-6).', badCnt);
        if strict, error(msg); else, warning(msg); end
    end
    
    fprintf('[VALIDATE] Basic DGP checks completed.\n');
end    
