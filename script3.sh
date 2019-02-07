for j in parameter_RMF012_ parameter_garnet parameter_RMF022_ parameter_RMF028_ parameter_RMF032_ parameter_gold2_
do
	for i in {0..199}
	do
		./main "$j$i" "18" "22" "$j Ar40" "5"
	done

	for i in {0..199}
	do
		./main "$j$i" "20" "16" "$j Ca36" "5"
	done

	for i in {0..199}
	do
		./main "$j$i" "16" "20" "$j S36" "5"
	done
done

