struct Hit
    motif::Int64 #index of motif
    strand::Int64 #0 or 1: palindromes??
    location::Int32
    score::Float32
end

function scan_seq(seq, motif, b_probs, hitsInSequences, seqnum, motnum, hit_thresh)
    tot_score::Float32=0
    m_max::Int64 = length(motif) #m_max is the num of motifs
    for m in 1:m_max
        pssm = deepcopy(motif[m])
        row_max = size(motif[m], 1)

        if(row_max > length(seq))
            continue
        end

        @fastmath @simd for r in 1:row_max
            @simd for c in 1:4
                @fastmath pssm[r, c] /= b_probs[c]
            end
        end

        #finally, scan the PSSM against the sequence:
        score::Float32 = 0
        posns = length(seq) - row_max + 1 #last possible beginning point for current motif in seq
        for n in 1:posns
            s=1
            for k in 1:row_max
                @fastmath s *= pssm[k, seq[n+k-1]+1]
            end
            score += s
            if(log(s) >= hit_thresh && seqnum != -1)
                push!(hitsInSequences[seqnum], Hit(motnum, m, n, s))
            end
        end
        @fastmath tot_score += (score/posns/m_max)
    end
   return tot_score
end
