for j in parameter_garnet_refined_ 
do
	for i in {0..199}
	do
		./main "$j$i" "6" "6" "$j C12" "3" "5"
	done

	
done

