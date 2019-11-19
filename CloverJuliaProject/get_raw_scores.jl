include("get_motifs.jl")
include("get_seqs.jl")
include("scan_seq.jl")
include("main.jl")
include("combine_score.jl")

singleStrandMotifs = Get_Single_Strand_Motifs("test.txt")
doubleStrandMotifs = Get_Double_Strand_Motifs(singleStrandMotifs, true)
#print("doubleStrandMotifs", doubleStrandMotifs)
sequenceNames = []
sequence = []
(sequence, sequenceNames) = get_seqs("fasta.txt", sequence, sequenceNames)
#print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
#print(length(sequence))
seq_info = seq_set_info(sequence)
base_probs = Vector{Float64}()
base_probs = get_base_probs(sequence, base_probs)

dp=Array{Float64}(undef,length(sequence[1])+1)

#p_value_temp = Array{Float64}(undef,1)
#seq_scores_temp = Array{Float64}(undef,1)
# temp_result = result(0,0.0,[],[])
results = []#Array{result,1}(undef,length(doubleStrandMotifs))
# print(methods(result))
motnum= 1 #temp
seqnum= 1 #temp
hit_thresh=6 #temp

seq_sco = []
raw_sco = []

print("\nlength of ds: ", length(doubleStrandMotifs))
for m in 1:length(doubleStrandMotifs)
    global raw_sco, temp_result, results
    # results = vcat(results,temp_result)
    # print(temp_result.motif_index)
    for s in 1:length(sequence)
        global seq_sco
        seq_sco = vcat(seq_sco, scan_seq(sequence[s], doubleStrandMotifs[m], base_probs, seqnum, motnum, hit_thresh))
    end
    raw_sco = vcat(raw_sco,Combine_Score(seq_sco[m]))
end

seq_sco = reshape(seq_sco,(length(doubleStrandMotifs),length(sequence)))
print("\nSequence score ", seq_sco)
print("\nRaw score ", raw_sco)

bg_info=[]
bg_seqs=[]


get_hits(sequence,doubleStrandMotifs,base_probs)
