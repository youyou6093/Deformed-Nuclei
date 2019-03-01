
for i in fsu020 iufsu parameter_nl3 tamufsua tamufsub tamufsuc
do
	./main "$i" "28" "22" "othermodels" "5"
	./main "$i" "22" "28" "othermodels" "5"
done


