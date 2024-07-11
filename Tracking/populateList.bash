#populating lookupList script (bozo)

#rm lingering files
cd /home/bhartsock/geantinoSimulator9000/build
rm data.csv
cd /home/bhartsock/lookupListV2
rm lookupList.csv
cd

for ((i=0; i<1000; i++)); do
rm nohup.out   
#Geant4
cd /home/bhartsock/geantinoSimulator9000/build
./LXe batch.mac
mv data.csv /home/bhartsock/lookupListV2
#C++
cd /home/bhartsock/lookupListV2
./lookupListV2
cd
done
