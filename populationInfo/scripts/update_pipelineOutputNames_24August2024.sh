  #!/usr/bin/env bash

### loop
  #sed '1d' /scratch/aob2x/DESTv2/populationInfo/scripts/update_DataFiles_26August2024_conversionTable.csv > \
  #/scratch/aob2x/update_DataFiles_26August2024_conversionTable.noHeader.csv


  while read p; do
    # p=$( cat /scratch/aob2x/update_DataFiles_26August2024_conversionTable.noHeader.csv | head -n1 )
    newName=$( echo ${p} | cut -d',' -f1 )
    oldName=$( echo ${p} | cut -d',' -f2 )
    echo $oldName
    wd=$( ls -d /project/berglandlab/DEST/dest_mapped/*/${oldName} )
    newwd=$( echo $wd | sed "s/${oldName}/$newName/g" )
    mkdir $newwd

    ls -d ${wd}/* | sed 's/ /\n/g' > /scratch/aob2x/tmp/${oldName}.fileList
    fileList=$( echo $fileList | tr ' ' '\n' )

    echo "24Aug2024 name change vocher" > $newwd/nameChangeVoucher_24Aug2024.txt


    while read -r line; do
        # line=$( echo ${fileList} | head -n1 )
        newFilename=$( echo ${line} | sed "s/${oldName}/${newName}/g" )
        echo "Renaming: $line   To: ${newFilename}" >> $wd/nameChangeVoucher_24Aug2024.txt
        mv ${line} ${newFilename}
    done < /scratch/aob2x/tmp/${oldName}.fileList

  done < /scratch/aob2x/update_DataFiles_26August2024_conversionTable.noHeader.csv
