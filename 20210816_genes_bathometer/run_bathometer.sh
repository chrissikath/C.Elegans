# quantify two loci that were overlapping with another gene in the annotation
# WBcel235 is reference genome
# loci for two genes from Daniel:
# I:7695681..7698222 gld1
# I:12483809..12486235 sygl1

for i in *.bed; do
    bathometer fpkm /home/udo/worms/*.idx :::  $i > quant_bathometer_all_genes_FPKM.tsv
done
