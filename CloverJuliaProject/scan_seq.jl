#struct hit
#    motif #index of motif
#    strand #0 or 1: palindromes??
#    location
#    score
#end

function scan_seq(seq, motif, b_probs, seqnum, motnum, hit_thresh)
    tot_score=0

    m_max=convert(Int64,size(motif)[1]/4)
    for m in 1:m_max
        pssm = zeros(Float64, (size(motif)[2], 4))#size(motif)[2] is row num
        row_max = convert(Int64,size(motif)[2])
        for row in 1:row_max
            for col in 1:4
                pssm[row,col] = motif[row+4*(m-1),col]
            end
        end

        if(size(pssm)[2] > length(seq))
            continue
        end

        for r in 1:row_max
            for c in 1:4
                pssm[r,c] /= b_probs[c]
            end
        end

        #finally, scan the PSSM against the sequence:
        score = 0
        posns = length(seq) - row_max +1 #last possible beginning point for current motif in seq

        for n in 1:posns
            s=1
            #for (uint k = 0; k < pssm.rows(); ++k)
               #s *= pssm[k][*(n+k)];
            #print(seq[1+1],"test\n")
            for k in 1:row_max-1
                s *= pssm[k,seq[n+k]+1]
            end
            score += s

            #hits[seqnum].push_back(hit(motnum, m, n - seq.begin(), s));
            #vector<vector<hit> > hits;  // motif locations in each sequence
            hits = Vector{}()
            if(seqnum!=0 && log(2,s)>= hit_thresh)
                push!(hits,hit(motnum,m,n-seq[1],s))#??
            end
        end

        tot_score += (score/posns/length(motif))

    end

    print(tot_score)
   return tot_score
end
