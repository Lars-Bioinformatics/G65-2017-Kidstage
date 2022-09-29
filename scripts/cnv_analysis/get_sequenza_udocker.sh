export PATH=/work/G65-2017-Kidstage/udocker:$PATH
udocker pull sequenza/sequenza
udocker create --name=sequenza sequenza/sequenza:latest