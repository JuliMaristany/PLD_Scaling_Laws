for f in mean*Q*profile*; do

        if [ -f "$f" ] # does file exist?

        then
			prot=$(echo $f |awk -F"_" '{print $2}')
			var=$(echo $f |awk -F"_" '{print $3}')
			echo ${prot}
			echo ${var}
			python3.9 denfitting.py ${prot} ${var} 12 12 0.001 0.2 0.3 0.9 
        fi
	#break
done
