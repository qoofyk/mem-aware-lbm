#!/usr/bin/env bash

# First created: 2020 Aug 17
# Last modified: 2020 Aug 17

# Author: Yuankun Fu
# email: qoofyk@gmail.com

root_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo $root_dir

if [[ $1 = crop ]]; then
    function crop {
        ./pdfcrop "$1" >/dev/null && mv "$(basename "$1" .pdf)"-crop.pdf "$1"
    }
else
    function crop { true; }
fi

# sequential exp
python3 ${root_dir}/seq_mflups.py -m "bridges" > bridges_seq_cube.csv
python3 ${root_dir}/render_mflups.py ./bridges_seq_cube.csv --xdata "dims" --xlabel "Side length of a 3D cube cavity"  --ylim "(5,16)" --no-xlog --no-ylog --legend "./legend-seq.csv"
${root_dir}/pdfcrop bridges_seq_cube.pdf bridges_seq_cube.pdf

# parallel cube inut exp
exptype="omp_cube"
for DIM in 112; do
  python3 ${root_dir}/omp_mflups.py -m "bridges" -d ${DIM} -e "${exptype}" > bridges_${exptype}_dim_${DIM}.csv
  python3 ${root_dir}/render_mflups.py ./bridges_${exptype}_dim_${DIM}.csv --ylim "(0,200)" --no-xlog --no-ylog
  ${root_dir}/pdfcrop bridges_${exptype}_dim_${DIM}.pdf bridges_${exptype}_dim_${DIM}.pdf
done

for DIM in 224 336 448 560 672 784 840; do
  python3 ${root_dir}/omp_mflups.py -m "bridges" -d ${DIM} -e ${exptype} > bridges_${exptype}_dim_${DIM}.csv
  python3 ${root_dir}/render_mflups.py ./bridges_${exptype}_dim_${DIM}.csv --ylim "(0,300)" --no-xlog --no-ylog
  ${root_dir}/pdfcrop bridges_${exptype}_dim_${DIM}.pdf bridges_${exptype}_dim_${DIM}.pdf
done

# fair cube, palabos change to -> step2's dimension
exptype="omp_fair_cube"
for DIM in 112; do
  python3 ${root_dir}/omp_mflups.py -m "bridges" -d ${DIM} -e "${exptype}" > bridges_${exptype}_dim_${DIM}.csv
  python3 ${root_dir}/render_mflups.py ./bridges_${exptype}_dim_${DIM}.csv --ylim "(0,200)" --no-xlog --no-ylog
  ${root_dir}/pdfcrop bridges_${exptype}_dim_${DIM}.pdf bridges_${exptype}_dim_${DIM}.pdf
done

for DIM in 224 336 448 560 672 784; do
  python3 ${root_dir}/omp_mflups.py -m "bridges" -d ${DIM} -e "${exptype}" > bridges_${exptype}_dim_${DIM}.csv
  python3 ${root_dir}/render_mflups.py ./bridges_${exptype}_dim_${DIM}.csv --ylim "(0,300)" --no-xlog --no-ylog
  ${root_dir}/pdfcrop bridges_${exptype}_dim_${DIM}.pdf bridges_${exptype}_dim_${DIM}.pdf
done