#!/bin/bash
HOST=bioftp.org
USER=FNP_KSJ
PASSWORD=FNp_71aofb
 
ftp -inv $HOST <<EOF
user $USER $PASSWORD

# cd /TBD160404_20160725/Sample_Ae_Chu_2
# mget * 
# cd /TBD160404_20160725/Sample_All1_Chu_3
# mget * 
# cd /TBD160404_20160726/Sample_Ae_Cho_2
# mget * 
# cd /TBD160404_20160726/Sample_Ae_Cho_3
# mget * 
# cd /TBD160404_20160726/Sample_Ae_Gae_2
# mget * 
# cd /TBD160404_20160726/Sample_Ae_Gae_3
# mget * 
# cd /TBD160404_20160726/Sample_Ae_Hu_2
# mget * 
# cd /TBD160404_20160726/Sample_Ae_Hu_3
# mget * 
# cd /TBD160404_20160726/Sample_Ae_Yeong_2
# mget * 
# cd /TBD160404_20160726/Sample_Ae_Yeong_3
# mget * 
# cd /TBD160404_20160726/Sample_All1_Cho_2
# mget * 
# cd /TBD160404_20160726/Sample_All1_Cho_3
# mget * 
# cd /TBD160404_20160726/Sample_All1_Chu_2
# mget * 
# cd /TBD160404_20160726/Sample_All1_Gae_2
# mget * 
# cd /TBD160404_20160726/Sample_All1_Gae_3
# mget * 
# cd /TBD160404_20160726/Sample_All1_Hu_2
# mget * 
# cd /TBD160404_20160726/Sample_All1_Hu_3
# mget * 
# cd /TBD160404_20160726/Sample_All1_Yeong2_160621
# mget * 
# cd /TBD160404_20160726/Sample_All1_Yeong3_160621
# mget * 
# cd /TBD160264_20160519/Sample_6
# get 6_2.fq.gz
# cd /TBD160264_20160519/Sample_8
# get 8_2.fq.gz 

# 08_03_2016
# cd /TBD160512_20160803/Sample_414F1-gea-1
# mget * 
# cd /TBD160512_20160803/Sample_414F1-gea-2
# mget * 
# cd /TBD160512_20160803/Sample_414F1-gea-3
# mget * 
# cd /TBD160512_20160803/Sample_414F1-gyeolcho-1
# mget * 
# cd /TBD160512_20160803/Sample_414F1-gyeolcho-2
# mget * 
# cd /TBD160512_20160803/Sample_414F1-gyeolcho-3
# mget * 
#  cd /TBD160512_20160803/Sample_414F1-gyeolhoo-1
# mget *  
# cd /TBD160512_20160803/Sample_414F1-gyeolhoo-2
# mget * 
# cd /TBD160512_20160803/Sample_414F1-gyeolhoo-3
# mget * 
# cd /TBD160512_20160803/Sample_414F1-young-1
# mget * 
# cd /TBD160512_20160803/Sample_414F1-young-2
# mget * 
# cd /TBD160512_20160803/Sample_414F1-young-3
# mget *
# cd /TBD160512_20160803/Sample_415F1-gea-1
# mget *
# cd /TBD160512_20160803/Sample_415F1-gea-2
# mget *
cd /TBD160512_20160803/Sample_415F1-gea-3
mget *
cd /TBD160512_20160803/Sample_415F1-gyeolcho-1
mget *
cd /TBD160512_20160803/Sample_415F1-gyeolcho-2
mget *
cd /TBD160512_20160803/Sample_415F1-gyeolcho-3
mget *
cd /TBD160512_20160803/Sample_415F1-gyeolhoo-1
mget *
cd /TBD160512_20160803/Sample_415F1-gyeolhoo-2
mget *
cd /TBD160512_20160803/Sample_415F1-gyeolhoo-3
mget *
cd /TBD160512_20160803/Sample_415F1-young-1
mget *
cd /TBD160512_20160803/Sample_415F1-young-2
mget *
cd /TBD160512_20160803/Sample_415F1-young-3
mget *


bye
EOF

