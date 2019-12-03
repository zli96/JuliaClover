include("globalVariable.jl")
using Random
rng= MersenneTwister();
Random.seed!(rng, randomSeed)
function get_fasta(fil)
    open(fil) do file
    titles=[];title="";seqs=Array{UInt8}[];seq=[]
        for s in eachline(file)
            if s[1]=='>'
                if title!=""
                    push!(titles,title)
                    push!(seqs,seq);seq=[]
                end
                title=SubString(s,2)
            elseif s[1]=='#'
                print("Comment ignored")
            else
                for x in s
                    if(DNA_to_number(x) != nothing)
                        append!(seq,DNA_to_number(x))
                    end
                end
            end
        end
        push!(titles,title)
        push!(seqs,seq)
    return seqs, titles
    end
end

function get_seqs(seq_file)
    myseqs=Array{UInt8}[];seq_names=[];n=""
    (myseqs,n)=get_fasta(seq_file)
    push!(seq_names,n)
    return myseqs, seq_names
    discard=pop!(myseqs)
end

function count_residues(seq, counts, ALPHSIZE)#error change seq[[]] to []
    if(counts == nothing || length(counts) < ALPHSIZE )
        counts = [0,0,0,0]
    end
    for x in 1:length(seq)
        if seq[x] < UInt(ALPHSIZE)
            counts[seq[x]+1] += 1
        end
    end
    return counts
end

function copy_masks(source, dest)
    for i in 1:length(source)
        if(source[i]==UInt8(4))
            dest[i] = UInt8(4)
        end
    end
    return dest
end

function get_base_probs(seq)
    probs = []
    counts=[]
    #count_residues(seq, counts, alphsize)
    temp_counts = count_residues(seq, counts, ALPHSIZE)
    counts  = temp_counts
    tot=0
    for x in counts
        tot=tot+x
    end
    #println("counts in get_base_probs: ", counts)
    for i = 1:ALPHSIZE
        push!(probs, counts[i]/tot)
    end
    return probs
end


mutable struct result
    motifIndex::Int64
    rawScore::Float32
    pValues::Array{Float32, 1}
    seqScores::Array{Float32, 1}
end

mutable struct seq_info
    num::Int64
    len::Int64
    gc::Float32
end
function init_seq_info(seqs)#seqs:[[]]
    num = length(seqs)
    len=0
    for s = 1:num
        len=len+length(seqs[s])
    end
    counts=[]
    for s = 1:length(seqs)
        temp_counts = count_residues(seqs[s], counts, ALPHSIZE)
        counts  = temp_counts
    end
    # println("in init_seq_info, counts: $counts")
    tot=0
    for x in counts
        tot=tot+x
    end
    gc = (counts[2]+counts[3])/tot
    return seq_info(num, len, gc)
end

function bg_fragment(bg_seqs, frag, len, frag_num, frag_tot)
    if(length(bg_seqs) == 0)
        return
    end
    b = 1 #select which bg seq
    r = rand(rng,1:frag_tot)
    for b in 1:length(bg_seqs)
        if(frag_num[b]>r)
            break
        end
        r -= frag_num[b]
    end

    #if(b == length(bg_seqs))
    #    return
    #end

    # p =  Vector{UInt}()
    posns = length(bg_seqs[b]) - len + 1
    p = bg_seqs[b][1]+rand(rng,1:posns)
    flag=true
    while (flag)
        p = bg_seqs[b][1]+rand(rng,1:posns-1)
        ind = 0
        # println(p," ",len)
        for i = p:p+len
            if bg_seqs[b][i] == ALPHSIZE
                ind=i
                break
            end
            if i==p+len
                ind=i
            end
        end
        if ind == p+len
            flag= false
        end
    end

    # push!(p,bg_seqs[b][1]+rand(1:posns))
    # for i in 1:length(p)+len
    #     if(bg_seqs[b][i]==UInt8(4)&& i!= length(p)+len)
    #         push!(p,bg_seqs[b][1]+rand(1:posns))
    #     end
    # end

    for i in p:p+len-1
        push!(frag,bg_seqs[b][i])
    end
    #println("frag: $frag")
    return
end

function is_significant(r)
    for p = 1:length(r.pValues)
        if (r.pValues[p] > PTHRESH && r.pValues[p] < 1-PTHRESH)
            return false
        end
    end
    return true
end



# bg_info = Array{}
# bg_files=[]
# for i = 1: length(bg_files)
#     bg_seqs = Array{Int32}[]
#     junk = Array{String}
#     get_seqs(bg_files[i], bg_seqs, junk)
#     bg_info=vcat(bg_info,init_seq_info(bg_seqs))
#     shuffle_bgseq(bg_seqs)
# end
# cts=[0,0,0,0]
# count_residues(s[1],cts,4)
