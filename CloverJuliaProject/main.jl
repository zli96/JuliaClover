function DNA_to_number(c)
    if c=='a' || c=='A'
        return UInt8(0)
    elseif c=='c' || c=='C'
        return UInt8(1)
    elseif c=='g' || c=='G'
        return UInt8(2)
    elseif c=='t' || c=='T'
        return UInt8(3)
    else
        return UInt8(4)
    end
end

function get_fasta(fil)
    open(fil) do file
    titles=[];title="";seqs=Array{UInt8}[];seq=[]
        for s in eachline(file)
            if s[1]=='>'
                if title!=""
                    titles=vcat(titles,title)
                    push!(seqs,seq);seq=[]
                end
                title=SubString(s,2)
            elseif s[1]=='#'
                print("Comment ignored")
            else
                for x in s
                    append!(seq,DNA_to_number(x))
                end
            end
        end
        titles=vcat(titles,title)
        push!(seqs,seq)
    return seqs, titles
    end
end

function get_seqs(seq_file, myseqs, seq_names)
    myseqs=Array{UInt8}[];seq_names=[];n=""
    (myseqs,n)=get_fasta(seq_file)
    seq_names=vcat(seq_names,n)
#     print(myseqs,seq_names)
    return myseqs, seq_names
    discard=pop!(myseqs)
end
# n=[];s=[]
#s,n=get_seqs("fasta.txt",s,n)
## print(out)
#print(s,n)
#if(isempty(out)==true)
#    print("Die : no sequences read\n")
#else
#    print("Set sequence info\n")
#end

alphsize=4
pthresh = 0.01
function count_residues(seq, counts, alphsize)
    if(length(counts)<alphsize)
        counts = [0,0,0,0]
    end
    for x = 1:length(seq)
        #print("A :",typeof(seq[x]))
        if seq[x] < UInt(alphsize)
            counts[seq[x]+1] = counts[seq[x]+1] + 1
        end
    end
    #print("updated",counts)
    return counts
end

function rand_test(myseqs, b_probs, motifs, losses)
    for m in 1:length(motifs)
        scores = Vector{Float64}()
        for s in 1:length(myseqs)
            a,b = scan_seq(myseqs[s], motifs[m], b_probs[s], 1, 1, 6)
            push!(scores,a)
        end
        raw = combine_score(scores)
        if (raw >= results[m].raw_score)#notes:no result[m]
          losses[m]+=1;
      end
    end
    return losses
end

function copy_masks(source, dest)
    for i in 1:length(source)
        if(source[i]==UInt8(4))
            dest[i] = UInt8(4)
        end
    end
end

function get_base_probs(seq, probs)
    counts=[]
    for s = 1:length(seq)
        #print("signature",seqs[s], counts, alphsize)
        counts=count_residues(seq[s], counts, alphsize)
    end
    #print("counts",counts)
    tot=0
    for x in counts
        tot=tot+x
    end
    for i = 1:alphsize
        probs=vcat(probs,counts[i]/tot)
    end
    return probs
end

#dp=Array{Float64}(undef,length(s[1])+1)  # INit dp line

mutable struct result
    motif_index::Int64
    raw_score::Float64
    pvalues::Array{Float64, 1}
    seq_scores::Array{Float64, 1}
end
function greaterThan(r1, r2)
    return r1.raw_score > r2.raw_score
end

mutable struct seq_info
    num::Int64
    len::Int64
    gc::Float64
end
function init_seq_info(seqs)
    num = length(seqs)
    len=0
    for s = 1:num
        len=len+length(seqs[s])
    end
    counts=[]
    print("\nseqs length ",length(seqs))
    for s = 1:length(seqs)
        #print("signature",seqs[s], counts, alphsize)
        counts=count_residues(seqs[s], counts, alphsize)
    end
    tot=0
    for x in counts
        tot=tot+x
    end
    gc = (counts[2]+counts[3])/tot
    return seq_info(num, len, gc)
end

function bg_fragment(bg_seqs, frag, len, frag_num, frag_tot)
    b = 1 #select which bg seq
    r = rand(1:frag_tot)
    for b in 1:length(bg_seqs)
        if(frag_num[b]>r)
            break
        end
        r -= frag_num[b]
    end

    if(b == length(bg_seqs))
        return
    end

    p =  Vector{UInt}()
    posns = length(bg_seqs[b]) - len + 1

    push!(p,bg_seqs[b][1]+rand(1:posns))
    for i in 1:length(p)+len
        if(bg_seqs[b][i]==UInt8(4)&& i!= length(p)+len)
            push!(p,bg_seqs[b][1]+rand(1:posns))
        end
    end

    for i in 1:length(p)
        push!(frag,length(p)+len)
    end
    return
end

function shuffle_bgseq(seqs,bg_seqs,ds_motifs,results)
    frag_nums = []
    for i = 1:length(seqs)
        t=[]

        for j = 1:length(bg_seqs)
            frags=0
            r=0
            for x in 1:length(bg_seqs[j])
                # global frag_nums
                if bg_seqs[j][x]==alphsize
                    r=0
                else
                    r=r+1
                    if r>length(seqs[i])
                        frags=frags+1
                    end
                end
                t=vcat(t,frags)
            end
            push!(frag_nums,t)
        end

    end
    # frag_nums = reshape(frag_nums,(length(seqs),length(bg_seqs)))
    print("\nfragnums",frag_nums)
    frag_tots = []
    for x in 1:length(frag_nums)
        tot=0
        for y in 1:length(frag_nums[x])
            tot=tot+frag_nums[x][y]
        end
        if tot==0
            print("\nDie function : Can't get fragments of control sequences to match all target sequences.")
        end
        push!(frag_tots,tot)
    end
    print("\nfrag tots",frag_tots)
    losses = zeros(length(ds_motifs))
    pvalues = []#Vector{Float64}(0,length(ds_motifs))
    shuffles = 1000
    for r = 1:shuffles
        r_seqs = []#Vector{Int64}
        b_probs = []#Vector{Float64}
        for s = 1:length(seqs)
            temp_seqs = []
            temp_probs = []
            bg_fragment(bg_seqs, temp_seqs, length(seqs[s]), frag_nums[s], frag_tots[s])
            push!(r_seqs,temp_seqs)
            copy_masks(seqs[s], r_seqs[s])
            get_base_probs(r_seqs[s], temp_probs)
            push!(b_probs,temp_probs)
        end
        losses = rand_test(r_seqs, b_probs, ds_motifs, losses)
    end

    for m = 1:length(ds_motifs)
        push!(results[m].pvalues,losses[m]/shuffles)
    end
    # return pvalues
end

# mutable struct hit
#     motif::Int64
#     strand::Int64
#     location::Int64
#     score::Float64
# end
function lessThan(h1, h2)
    return h1.location < h2.location
end
function is_significant(r)
    for p = 1:length(r.pvalues)
        if (r.pvalues[p] > pthresh) && (r.pvalues[p] < 1-pthresh)
            return false
        end
    end
    return true
end

function get_hits(seqs,ds_motifs,base_probs,results,hits)
    resize!(hits,length(seqs))
    for m = 1:length(ds_motifs)
        if is_significant(results[m])
            for s = 1:length(seqs)
                scan_seq(seqs[s], ds_motifs[m], base_probs, s, m, 6)
            end
        end
    end
end

# bg_info = Array{}
bg_files=[]
for i = 1: length(bg_files)
    bg_seqs = Array{Int32}[]
    junk = Array{String}
    get_seqs(bg_files[i], bg_seqs, junk)
    bg_info=vcat(bg_info,init_seq_info(bg_seqs))
    shuffle_bgseq(bg_seqs)
end
# cts=[0,0,0,0]
# count_residues(s[1],cts,4)
