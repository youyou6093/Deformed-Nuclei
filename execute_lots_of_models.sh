for j in parameter_garnet parameter_RMF012_ parameter_RMF022_ parameter_RMF028_ parameter_RMF032_
do
	for i in {0..199}
	do
		./main "$j$i" "54" "78" "$j Xe132" "7" "5"
	done

	
done

