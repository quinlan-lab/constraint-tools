# usage ./lift.sh input input.version output.version "FLAGS"
#
# e.g: lift.sh Homo_sapiens.GRCh37.61.gtf hg19 hg18 -gff
# will map hg19 features to hg18 coordinates.
if [ -z $3 ]
then
    echo "usage: lift Homo_sapiens.GRCh37.61.gtf hg19 hg18 -gff"
    exit 1
fi

DIR=~/.lift/

mkdir -p $DIR

INPUT=$1
FROM_VERSION=$2
TO_VERSION=$3
FLAGS=$4

TO_CAPS=`echo $TO_VERSION | sed 's/\<./\u&/'`

if [[ ! -f $DIR/liftOver ]]
then
    wget -O $DIR/liftOver http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
    chmod +x $DIR/liftOver
fi

CHAINF=${FROM_VERSION}To${TO_CAPS}.over.chain
if [[ ! -f $DIR/$CHAINF ]]
then
   wget -O $DIR/$CHAINF.gz http://hgdownload.cse.ucsc.edu/goldenPath/${FROM_VERSION}/liftOver/${CHAINF}.gz
   gunzip $DIR/${CHAINF}.gz
fi

OUTPUT=$INPUT.$TO_VERSION

$DIR/liftOver $FLAGS -minMatch=0.8 $INPUT $DIR/$CHAINF $OUTPUT $OUTPUT.unmapped

wc -l $OUTPUT
wc -l $OUTPUT.unmapped