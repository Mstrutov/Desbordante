git lfs pull
if [ $? -eq 0 ]; then
    echo "PULLED VIA GIT LFS"
else
    echo "FAILED TO PULL VIA GIT LFS, RETRIEVING FROM 'Desbordante-Data' REPOSITORY"
    git clone https://github.com/Mstrutov/Desbordante-Data.git
    mv Desbordante-Data/datasets.zip datasets/datasets.zip
    rm -rf Desbordante-Data
fi
