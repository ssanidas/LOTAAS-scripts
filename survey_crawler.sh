#/bin/bash

#========================#
#GBT driftscan survey 350#
#========================#

#wget the website
wget -O - http://astro.phys.wvu.edu/GBTdrift350 >GBTds350web.txt

#get the source entries
cat GBTds350web.txt |grep "^<td><a href" >GBTds350refined.txt

#psr names
cat GBTds350refined.txt |awk -F">" '{print $3}'|awk -F"<" '{print $1}' >GBTds350names.txt

#psr periods in ms
cat GBTds350refined.txt |awk -F"> " '{print $2}'|awk -F" <" '{print $1}' >GBTds350periods.txt

#psr DMs
cat GBTds350refined.txt |awk -F"> " '{print $3}'|awk -F" <" '{print $1}' >GBTds350DM.txt

#psr RA
cat GBTds350names.txt|awk '{print substr($0,2,2)":"substr($0,4,2)}' >GBTds350RA.txt

#psr DEC
cat GBTds350names.txt|awk '{if (length($0) == 8) print substr($0,6,3);else print substr($0,6,3)":"substr($0,9,2)}' >GBTds350DEC.txt

#creating the survey tag
nsource=`cat GBTds350names.txt |wc -l`
for i in `seq 1 $nsource`;do
    echo "GBTds350" >> GBTds350tag.txt
    done

#creating the combined file
paste GBTds350names.txt GBTds350RA.txt GBTds350DEC.txt GBTds350periods.txt GBTds350DM.txt GBTds350tag.txt>GBTds350parsed.txt 

#cleaning up
rm GBTds350web.txt GBTds350refined.txt GBTds350names.txt GBTds350periods.txt GBTds350DM.txt GBTds350RA.txt GBTds350DEC.txt GBTds350tag.txt

#============#
#PALFA survey#
#============#

#wget the website
wget -O - http://www.naic.edu/~palfa/newpulsars/ >PALFAweb.txt

#remove html comments
cat PALFAweb.txt|sed -e :a -re 's/<!--.*?-->//g;/<!--/N;//ba'  >PALFAwebclean.txt

#remove empty spaces
cat PALFAwebclean.txt|tr -d " \t\n\r" >PALFAwebstrippedclean.txt

#get the number of source entries
nsource=`cat PALFAwebclean.txt |grep "<tr align=center>"|wc -l`
nsource=$(($nsource+1))

#looping through the entries
for i in `seq 3 $nsource`;do
    psrname=`cat PALFAwebstrippedclean.txt|awk -F"<tralign=center>" '{print $'$i'}'|awk -F"<td>" '{print $3}'|awk -F"</td>" '{print $1}'`
    if [[ "$psrname" == *"href"* ]];then
	echo  "$psrname"|grep -Po "(?<=\>).*?(?=\<)" >>PALFAnames.txt
    else
	echo "$psrname"  >>PALFAnames.txt #get the psr name (modified for the case where the psr name is a link)
    fi
    cat PALFAwebstrippedclean.txt|awk -F"<tralign=center>" '{print $'$i'}'|awk -F"<td>" '{print $4}'|awk -F"</td>" '{print $1}'|sed 's/<b>//g'|sed 's\</b>\\g'|sed 's/~//g' >>PALFAperiod.txt #get the period in ms
    cat PALFAwebstrippedclean.txt|awk -F"<tralign=center>" '{print $'$i'}'|awk -F"<td>" '{print $5}'|awk -F"</td>" '{print $1}'|sed 's/<b>//g'|sed 's\</b>\\g' >>PALFADM.txt #get the DM
    done

#psr RA
cat PALFAnames.txt|awk '{print substr($0,2,2)":"substr($0,4,2)}' >PALFARA.txt

#psr DEC
cat PALFAnames.txt|awk '{if (length($0) == 8) print substr($0,6,3);else print substr($0,6,3)":"substr($0,9,2)}' >PALFADEC.txt

# merging the files
paste PALFAnames.txt PALFARA.txt PALFADEC.txt PALFAperiod.txt PALFADM.txt >PALFAparsed.txt 


#cleaning up
rm PALFAwebclean.txt PALFAnames.txt PALFAperiod.txt PALFADM.txt PALFAwebstrippedclean.txt PALFAweb.txt PALFARA.txt PALFADEC.txt


#======================#
#AO 327MHz drift survey#
#======================#

#wget the website
wget -O - http://www.naic.edu/~deneva/drift-search/ >AO327web.txt

#remove html comments
cat AO327web.txt|sed -e :a -re 's/<!--.*?-->//g;/<!--/N;//ba'  >AO327webclean.txt

#remove empty spaces
cat AO327webclean.txt|tr -d " \t\n\r" >AO327webstrippedclean.txt

#get the number of source entries
nsource=`cat AO327webclean.txt |grep "<tr align=center>"|wc -l`
nsource=$(($nsource+1))

#looping through the entries
for i in `seq 3 $nsource`;do
    psrname=`cat AO327webstrippedclean.txt|awk -F"<tralign=center>" '{print $'$i'}'|awk -F"<td>" '{print $3}'|awk -F"</td>" '{print $1}'`
    if [[ "$psrname" == *"href"* ]];then
        echo  "$psrname"|grep -Po "(?<=\>).*?(?=\<)"|sed 's/*//g' >>AO327names.txt
    else
        echo "$psrname"|sed 's/*//g'  >>AO327names.txt #get the psr name (modified for the case where the psr name is a link)
    fi
    psrperiod=`cat AO327webstrippedclean.txt|awk -F"<tralign=center>" '{print $'$i'}'|awk -F"<td>" '{print $4}'|awk -F"</td>" '{print $1}'|sed 's/'"<b>"'//g'|sed 's/'"<\/b>"'//g'|sed 's/'"~"'//g'`
    if [[ "$psrperiod" == "&nbsp" ]];then
	echo ""  >>AO327period.txt
    else
	echo "$psrperiod" >>AO327period.txt
    fi  #get the period in ms    
    cat AO327webstrippedclean.txt|awk -F"<tralign=center>" '{print $'$i'}'|awk -F"<td>" '{print $5}'|awk -F"</td>" '{print $1}'|sed 's/<b>//g'|sed 's\</b>\\g' >>AO327DM.txt #get the DM
done

#psr RA
cat AO327names.txt|awk '{print substr($0,2,2)":"substr($0,4,2)}' >AO327RA.txt

#psr DEC
cat AO327names.txt|awk '{if (length($0) == 8) print substr($0,6,3);else print substr($0,6,3)":"substr($0,9,2)}' >AO327DEC.txt

# merging the files 
paste AO327names.txt AO327RA.txt AO327DEC.txt AO327period.txt AO327DM.txt >AO327parsed.txt

#cleaning up
rm AO327webclean.txt AO327names.txt AO327period.txt AO327DM.txt AO327webstrippedclean.txt AO327web.txt AO327RA.txt AO327DEC.txt


#==================#
#ATNF PSR Catalogue#
#==================#

#IMPORTANT!
#the website will need adjustment when the catalogue will be updated

#wget the catalogue
wget -O - "http://www.atnf.csiro.au/research/pulsar/psrcat/proc_form.php?version=1.53&JName=JName&RaJ=RaJ&DecJ=DecJ&P0=P0&DM=DM&startUserDefined=true&c1_val=&c2_val=&c3_val=&c4_val=&sort_attr=jname&sort_order=asc&condition=&pulsar_names=&ephemeris=short&coords_unit=raj%2Fdecj&radius=&coords_1=&coords_2=&style=Short+csv+without+errors&no_value=*&fsize=3&x_axis=&x_scale=linear&y_axis=&y_scale=linear&state=query&table_bottom.x=73&table_bottom.y=13" >ATNFweb.txt

#pick up only the source entries and remove any empty lines
sed '/^[^0-9]/d' ATNFweb.txt |sed '/^$/d' >ATNFclean.txt

#get the psr names
cat ATNFclean.txt|awk -F";" '{print $2}' >ATNFnames.txt

#get the RA
cat ATNFclean.txt|awk -F";" '{print $3}' >ATNFRA.txt

#get the DEC
cat ATNFclean.txt|awk -F";" '{print $4}' >ATNFDEC.txt

#get the period
cat ATNFclean.txt|awk -F";" '{print $5}' >ATNFperiod.txt

#get the DM
cat ATNFclean.txt|awk -F";" '{print $6}'|sed 's/*//g' >ATNFDM.txt

#merge the files
paste ATNFnames.txt ATNFRA.txt ATNFDEC.txt ATNFperiod.txt ATNFDM.txt > ATNFparsed.txt

#cleaning up
rm ATNFweb.txt ATNFclean.txt ATNFDM.txt ATNFperiod.txt ATNFDEC.txt ATNFRA.txt ATNFnames.txt

