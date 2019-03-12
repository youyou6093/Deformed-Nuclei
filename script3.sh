for j in parameter_RMF012_ parameter_garnet parameter_RMF022_ parameter_RMF028_ parameter_RMF032_ parameter_gold2_
do
	for i in {0..199}
	do
		./main "$j$i" "20" "28" "$j Ca48" "5"
	done

	
done

