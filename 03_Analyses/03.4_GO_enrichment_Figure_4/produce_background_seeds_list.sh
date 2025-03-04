# Produce list of seeds for background (Gene Universe) -- enrichment analysis

# V. cardui
cd cardui
grep -v '#' cardui.emapper.annotations | cut -f2,10 | sort -u > cardui.seed-GO.list
bash scripts/enrichment/genetoGO.sh cardui.seed-GO.list cardui.seed-GO.2cols.list
awk '{print $2"\t"$1}' cardui.seed-GO.2cols.list > cardui.GO-seed.2cols.list
rm cardui.seed-GO.2cols.list
cut -f1 cardui.seed-GO.list > cardui.seeds.list


# V. atalanta
cd atalanta
grep -v '#' atalanta.emapper.annotations | cut -f2,10 | sort -u > atalanta.seed-GO.list
bash scripts/enrichment/genetoGO.sh atalanta.seed-GO.list atalanta.seed-GO.2cols.list
awk '{print $2"\t"$1}' atalanta.seed-GO.2cols.list > atalanta.GO-seed.2cols.list
rm atalanta.seed-GO.2cols.list
cut -f1 atalanta.seed-GO.list > atalanta.seeds.list


# A. io
cd aglaisio
grep -v '#' aglaisio.emapper.annotations | cut -f2,10 | sort -u > aglaisio.seed-GO.list
bash scripts/enrichment/genetoGO.sh aglaisio.seed-GO.list aglaisio.seed-GO.2cols.list
awk '{print $2"\t"$1}' aglaisio.seed-GO.2cols.list > aglaisio.GO-seed.2cols.list
rm aglaisio.seed-GO.2cols.list
cut -f1 aglaisio.seed-GO.list > aglaisio.seeds.list


# V. tameamea
cd tameamea
grep -v '#' tameamea.emapper.annotations | cut -f2,10 | sort -u > tameamea.seed-GO.list
bash scripts/enrichment/genetoGO.sh tameamea.seed-GO.list tameamea.seed-GO.2cols.list
awk '{print $2"\t"$1}' tameamea.seed-GO.2cols.list > tameamea.GO-seed.2cols.list
rm tameamea.seed-GO.2cols.list
cut -f1 tameamea.seed-GO.list > tameamea.seeds.list