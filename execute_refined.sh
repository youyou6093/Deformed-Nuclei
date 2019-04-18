for j in parameter_RMF028_refined_
do
	for i in {0..199}
	do
		./main "$j$i" "20" "28" "$j Ca48" "5"
		./main "$j$i" "6" "6" "$j C12" "5"
	done

	
done

