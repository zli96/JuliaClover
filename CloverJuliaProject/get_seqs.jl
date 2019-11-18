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

#function get_seqs(seq_file, myseqs, seq_names)
#    myseqs=Array{UInt8}[]
#    seq_names=[]
#    n=""
#    (myseqs,n)=get_fasta(seq_file)
#    seq_names=vcat(seq_names,n)
#    print(myseqs,seq_names)
#    discard=pop!(myseqs)
#end

function get_base_probs(seq, probs)
    counts=[]
    count_residues(seq, counts, alphsize)
    tot=0
    for x in counts
        tot=tot+x
    end
    for i = 1:alphsize
        probs=vcat(probs,counts[i]/tot)
    end
end
