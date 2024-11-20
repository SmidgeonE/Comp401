if [ $# -eq 0 ]; then
    echo "Usage: $0 run -> This submits all jobs in the bc4 directory"
    echo "Usage: $0 check -> This checks the status of all jobs submitted by oz21652"
    echo "Usage: $0 cancel -> This cancels all jobs submitted by oz21652"

    exit 0
fi


case $1 in
    run)
        for file in /user/home/oz21652/Comp401/bc4/*.sh; do
            sbatch "$file"
        done

        ;;
    check)
        squeue -u oz21652
        ;;
    cancel)
        scancel -u oz21652
        ;;
esac