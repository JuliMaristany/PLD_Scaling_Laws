for f in *.profile; do

        if [ -f "$f" ] # does file exist?

        then
  			echo $f
			echo "tail_${f}"
			tail -n 5050 $f > "tail_${f}"
			split -l 1010 "tail_${f}" "tail_${f%"profile"}" 
        fi
	#break
done
