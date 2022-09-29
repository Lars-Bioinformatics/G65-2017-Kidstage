export PATH=/work/G65-2017-Kidstage/udocker:$PATH

udocker pull quay.io/wtsicgp/dockstore-cgpmap:3.1.4
udocker create --name=cgpmap quay.io/wtsicgp/dockstore-cgpmap:3.1.4

udocker pull quay.io/wtsicgp/dockstore-cgpwgs:2.1.1
udocker create --name=cgpwgs quay.io/wtsicgp/dockstore-cgpwgs:2.1.1

udocker pull quay.io/wtsicgp/cgpbattenberg:3.7.1
udocker create --name=cgpbattenberg quay.io/wtsicgp/cgpbattenberg:3.7.1

# udocker pull godlovedc/cgpbattenberg:3.5.3
# udocker create --name=cgpbattenberg353 godlovedc/cgpbattenberg:3.5.3
udocker pull quay.io/wtsicgp/cgpbattenberg:3.5.3
udocker create --name=cgpbattenberg353 quay.io/wtsicgp/cgpbattenberg:3.5.3