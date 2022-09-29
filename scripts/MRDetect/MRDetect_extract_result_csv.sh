# for dir in G65*; do
#     cd $dir
#     mkdir -p ${dir%_ni*vs*}_results
#     cp *csv ${dir%_ni*vs*}_results
#     cd ..
# done

# zip -r ichorCNA_results.zip */*csv
zip -r MRDetect_results.zip */*csv