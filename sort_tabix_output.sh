#!/usr/bin/bash

temp=$1
out=$2
tabix_str=${3:--b2 -e2 -S1}
do_sort=${4:-1}
do_tabix=${5:-1}


if [[ "$do_sort"x == "1"x ]]; then
	echo 'Sorting output file.'
	head -n1 ${temp} > ${out}

	tail -n+2 ${temp} | \
	sort -nk1 -nk2 >> ${out} && \
	rm ${temp} || \
	( tail -n+2 ${temp} | \
	sort -T ${out%/*} -nk1 -nk2 >> ${out} && \
	rm ${temp} ) && \
	echo 'Sort successful.' || \
	( echo 'Sort failed.'; exit 1 )
	
fi


if [[ "$do_tabix"x == "1"x ]]; then
	echo 'Bgzip/tabix-ing.'
	bgzip -f ${out} && \
	tabix -f ${tabix_str} ${out}.gz
fi




# if [ ! -f "${temp}" ] ; then
# 	echo 'Sort successful.'

# 	if [[ "$do_tabix"x == "1"x ]]; then
# 		echo 'Bgzip/tabix-ing.'
# 		bgzip -f ${out} && \
# 		tabix -f ${tabix_str} ${out}.gz
# 	fi

# else
# 	echo 'Sort failed.'

# 	if [[ "$do_tabix"x == "1"x ]]; then
# 		echo 'Unable to bgzip/tabix output file.'
# 	fi	

# fi


# echo 'Sorting output file.'
# head -n1 ${temp} > ${out}

# tail -n+2 ${temp} | \
# sort -nk1 -nk2 >> ${out} && \
# rm ${temp} || \
# ( tail -n+2 ${temp} | \
# sort -T ${out%/*} -nk1 -nk2 >> ${out} && \
# rm ${temp} )

# if [ ! -f "${temp}" ] ; then
# 	echo 'Sort successful. Bgzip/tabix-ing.'

# 	bgzip -f ${out} && \
# 	tabix -f ${tabix_str} ${out}.gz

# else
# 	echo 'Sort failed; Unable to bgzip/tabix output file.'
# fi

