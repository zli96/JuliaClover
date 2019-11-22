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

function get_seqs(seq_file)
    myseqs=Array{UInt8}[];seq_names=[];n=""
    (myseqs,n)=get_fasta(seq_file)
    seq_names=vcat(seq_names,n)
    return myseqs, seq_names
    discard=pop!(myseqs)
end

alphsize=4
pthresh = 0.01
function count_residues(seq, counts, alphsize)
    if(length(counts) < alphsize)
        counts = [0,0,0,0]
    end
    for x = 1:length(seq)
        if seq[x] < UInt(alphsize)
            counts[seq[x]+1] = counts[seq[x]+1] + 1
        end
    end
    #print("counts: ", counts)
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

function get_base_probs(seq, probs)
    counts=[]
    for s = 1:length(seq)
        #print("get_base_probs: ", seq[s], counts, alphsize)
        counts=count_residues(seq[s], counts, alphsize)
    end
    tot=0
    for x in counts
        tot=tot+x
    end
    #println("counts: ", counts)
    #println("tot: ", tot)
    for i = 1:alphsize
        push!(probs, counts[i]/tot)
    end
    return probs
end

#dp=Array{Float64}(undef,length(s[1])+1)  # INit dp line

mutable struct result
    motifIndex::Int64
    rawScore::Float64
    pValues::Array{Float64, 1}
    seqScores::Array{Float64, 1}
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
    for s = 1:length(seqs)
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



function lessThan(h1, h2)
    return h1.location < h2.location
end
function is_significant(r)
    for p = 1:length(r.pValues)
        if (r.pValues[p] > pthresh) && (r.pValues[p] < 1-pthresh)
            return false
        end
    end
    return true
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
