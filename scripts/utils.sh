function pair_epiread {
  sort -k2,2 -k3,3n "$1" |
      awk 'BEGIN{qname="";rec=""}
         qname==$2{print rec"\t"$5"\t"$6"\t"$7"\t"$8;qname=""}
         qname!=$2{qname=$2;rec=$1"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8;pair=$3}' >"$2"
}

function pair_epiread2 {
  sort -k2,2 -k3,3n "$1" |
      awk 'BEGIN{qname="";rec=""}
         qname==$2{print rec"\t"$5"\t"$6;qname=""}
         qname!=$2{qname=$2;rec=$1"\t"$4"\t"$5"\t"$6;pair=$3}' >"$2"
}
