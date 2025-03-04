# Run eggNOG-mapper to retrieve protein ortholog hits

# V. cardui
emapper.py -i V_card2.1_protein.faa -o cardui --cpu 20 --itype proteins

# V. atalanta
emapper.py -i V_atal_protein.faa -o atalanta --cpu 20 --itype proteins

# A. io
emapper.py -i A_io_protein.faa -o aglaisio --cpu 20 --itype proteins

# V. tameamea
emapper.py -i V_tameamea_protein.faa -o tameamea --cpu 20 --itype proteins
