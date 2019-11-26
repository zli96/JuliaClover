include("get_motifs.jl")
include("get_seqs.jl")
include("scan_seq.jl")
include("main.jl")
include("combine_score.jl")
using DataFrames
using CSV

function print_result(sequenceFileName, motifFileName, seq_info, motifSize, results)
    println("Sequence file: $sequenceFileName ($(seq_info.num) sequences, $(seq_info.len) bp, $(seq_info.gc) GC content)")
    println("Motif file: $motifFileName ($motifSize motifs)")

    for i in 1:length(results)
        results[i].motifIndex = i
    end
    sort!(results, by = x -> x.rawScore)

    println("\n ************* Over and under represented motifs *************")
    summary = DataFrame(Motif = String[], RawScore = Float64[], PValue = String[])
    for i in 1:length(results)
        if(is_significant(results[i]))
            motifName = motifNames[i]
            motifRawScore = results[i].rawScore
            pValueAsString = ""
            for j in results[i].pValues
                pValueAsString = string(pValueAsString, string(j))
            end
            push!(summary, [motifName, motifRawScore, pValueAsString])
        end

    end
    display(summary)
    CSV.write("summary.csv", summary)

    println("*** Motif Instances with Score >= 6:")
    for i in 1:length(hitsInSequences)
        hitFrame = DataFrame(Motif = String[], Location = String[], Strand = String[], Sequence = String[], Score = Float64[])
        for j in 1:length(hitsInSequences[i])
            hit = hitsInSequences[i][j]
            strand = ""
            if(hit.strand == 1)
                strand = "+"
            elseif(hit.strand == 2)
                strand = "-"
            else
                print("Something wrong here! hit#: $j, sequence#: $i")
            end

            motifWidth = length(singleStrandMotifs[hit.motif])
            location = String("$(hit.location) - $(hit.location + motifWidth - 1)")
            motifString = ""
            for w in 0:(motifWidth-1)
                motifString = string(motifString, number_to_DNA(sequence[i][(hit.location + w)]))
            end
            push!(hitFrame, [motifNames[hitsInSequences[i][j].motif], location, strand, motifString, log(hit.score)])
        end
        display(hitFrame)
        fileName = string("hitInSequence", i)
        fileName = string(fileName, ".csv")
        CSV.write(fileName, hitFrame)
    end
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
                push!(t,frags)
            end
            push!(frag_nums,t)
        end
    end
    println("computed frag_num")
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
    println("computed frag_tots")
    losses = zeros(length(ds_motifs))
    pValues = []#Vector{Float64}(0,length(ds_motifs))
    shuffles = 2
    for r = 1:shuffles
        println("The $r th shuffle")
        r_seqs = []#Vector{Int64}
        b_probs = []#Vector{Float64}
        for s = 1:length(seqs)
            temp_seqs = []
            bg_fragment(bg_seqs, temp_seqs, length(seqs[s]), frag_nums[s], frag_tots[s])
            push!(r_seqs,temp_seqs)
            r_seqs[s] = copy_masks(seqs[s], r_seqs[s])
            temp_probs = get_base_probs(r_seqs[s])
            push!(b_probs, temp_probs)
        end
        losses = rand_test(r_seqs, b_probs, ds_motifs, losses)
    end

    for m = 1:length(ds_motifs)
        push!(results[m].pValues,losses[m]/shuffles)
    end
    # return pValues
end

function rand_test(myseqs, b_probs, motifs, losses)
    for m in 1:length(motifs)
        scores = Vector{Float64}()
        for s in 1:length(myseqs)
            a = scan_seq(myseqs[s], motifs[m], b_probs[s], hitsInSequences, s, m, 6)
            push!(scores,a)
        end
        raw = combine_score(scores)
        if (raw >= results[m].rawScore)#notes:no result[m]
          losses[m]+=1;
      end
    end
    return losses
end

function get_hits(seqs,ds_motifs,base_probs,results,hits)
    resize!(hits,length(seqs))
    for m in 1:length(ds_motifs)
        if (is_significant(results[m]))
            for s in 1:length(seqs)
                scan_seq(seqs[s], ds_motifs[m], base_probs[s], hitsInSequences, s, m, 6)
            end
        end
    end
end


bg_info=[]
bg_seqs=[]
bg_sequenceNames=[]
(bg_seqs, junk) = get_seqs("hs_chr20 - Copy.mfa")
sequenceFileName = "hs_dopamine_nr.fa - Copy.txt"
motifFileName = "jaspar2009core - Copy.txt"

#(bg_seqs, junk) = get_seqs("bg_seqs.txt")
#sequenceFileName = "fasta.txt"
#motifFileName = "test.txt"


motnum= 1 #temp
hit_thresh=6 #temp

(singleStrandMotifs, motifNames) = Get_Single_Strand_Motifs(motifFileName)
doubleStrandMotifs = Get_Double_Strand_Motifs(singleStrandMotifs, true)
#println(doubleStrandMotifs)
sequenceNames = []
sequence = []
(sequence, sequenceNames) = get_seqs(sequenceFileName)
#println(sequence[1])
seqInfo = init_seq_info(sequence)
base_probs = []
for s in 1:length(sequence)
    temp_base_probs = get_base_probs(sequence[s])
    push!(base_probs, temp_base_probs)
end

dp=Array{Float64}(undef,length(sequence[1])+1)

results = Vector{result}(undef,length(doubleStrandMotifs))

hitsInSequences = [Vector{Hit}() for i in 1:length(sequence)]
for m in 1:length(doubleStrandMotifs)
    sequenceScores = []
    #println("scanning motif $m")
    for s in 1:length(sequence)
        a = scan_seq(sequence[s], doubleStrandMotifs[m], base_probs[s], hitsInSequences, s, motnum, hit_thresh)
        push!(sequenceScores, a)
    end
    rawScore = combine_score(sequenceScores)
    results[m] = result(0, rawScore, [], sequenceScores)
end

bg_info = vcat(bg_info,init_seq_info(bg_seqs))
println("starting to shuffle!")
@time p_vals = shuffle_bgseq(sequence, bg_seqs, doubleStrandMotifs, results)

get_hits(sequence, doubleStrandMotifs, base_probs, results, hitsInSequences)

print_result(sequenceFileName, motifFileName, seqInfo, length(doubleStrandMotifs), results)
