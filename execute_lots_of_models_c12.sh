for j in parameter_garnet_refined_ parameter_rmf022_refined_ parameter_rmf028_refined_ parameter_rmf032_refined_ 
do
	for i in {0..199}
	do
		./main "$j$i" "6" "6" "$j C12" "3"
	done

	
done

