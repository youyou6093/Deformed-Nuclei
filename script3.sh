for j in parameter_RMF012_ parameter_garnet parameter_RMF022_ parameter_RMF028_ parameter_RMF032_ parameter_gold2_
do
	for i in {0..199}
	do
		./main "$j$i" "28" "22" "$j Ni50" "5"
	done

	for i in {0..199}
	do
		./main "$j$i" "22" "28" "$j Ti50" "5"
	done
done
