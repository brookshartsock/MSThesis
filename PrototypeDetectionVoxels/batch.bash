#written with the help of GPT4
input_directories=(
  "/home/bhartsock/Desktop/pTube_data/FD/Na-22"
  "/home/bhartsock/Desktop/pTube_data/FD/Cs-137"
  "/home/bhartsock/Desktop/pTube_data/FD/Mn-54"
  "/home/bhartsock/Desktop/pTube_data/BL1/Na-22"
  "/home/bhartsock/Desktop/pTube_data/BL1/Cs-137"
  "/home/bhartsock/Desktop/pTube_data/BL1/Mn-54"
  "/home/bhartsock/Desktop/pTube_data/BL2/Na-22"
  "/home/bhartsock/Desktop/pTube_data/BL2/Cs-137"
  "/home/bhartsock/Desktop/pTube_data/BL2/Mn-54"
  "/home/bhartsock/Desktop/pTube_data/HL/Na-22"
  "/home/bhartsock/Desktop/pTube_data/HL/Cs-137"
  "/home/bhartsock/Desktop/pTube_data/HL/Mn-54"
)

output_files=(
  "FDNa_peaksSiPM.csv"
  "FDCs_peaksSiPM.csv"
  "FDMn_peaksSiPM.csv"
  "BL1Na_peaksSiPM.csv"
  "BL1Cs_peaksSiPM.csv"
  "BL1Mn_peaksSiPM.csv"
  "BL2Na_peaksSiPM.csv"
  "BL2Cs_peaksSiPM.csv"
  "BL2Mn_peaksSiPM.csv"
  "HLNa_peaksSiPM.csv"
  "HLCs_peaksSiPM.csv"
  "HLMn_peaksSiPM.csv"
)

for ((i=0; i<${#input_directories[@]}; i++)); do
  echo "Input (directory):${input_directories[i]}" > peakIn.txt
  echo "Output (file):${output_files[i]}" >> peakIn.txt
  echo "Running ./landauFitPeaks with peakIn.txt..."
  ./landauFitPeaks peakIn.txt
done
