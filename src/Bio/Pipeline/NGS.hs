{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE GADTs             #-}
{-# LANGUAGE OverloadedLists   #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}

module Bio.Pipeline.NGS
    ( BWAOpts
    , BWAOptSetter
    , bwaCores
    , bwaSeedLen
    , bwaMaxMis
    , bwaReadTrim
    , defaultBWAOpts
    , bwaMkIndex
    , bwaAlign
    , filterBam
    , removeDuplicates
    , bam2Bed
    , sortedBam2BedPE
    , mergeReplicatesBed

    , STAROpts
    , STAROptSetter
    , starCmd
    , starCores
    , starSort
    , starTmpDir
    , starMkIndex
    , starAlign

    , RSEMOpts
    , RSEMOptSetter
    , rsemPath
    , rsemCores
    , rsemSeed
    , rsemMkIndex
    , rsemQuant
    ) where

import           Bio.Data.Bam              (bamToBed, readBam, runBam,
                                            sortedBamToBedPE)
import           Bio.Data.Bed              (BED3 (..), BEDLike (..), toLine)
import           Bio.Data.Experiment.Types
import           Conduit
import           Control.Lens
import           Control.Monad             (forM)
import           Control.Monad.State.Lazy
import           Data.Conduit.Zlib         (gzip)
import           Data.Int                  (Int32)
import           Data.Maybe                (fromJust)
import qualified Data.Text                 as T
import           Shelly                    (escaping, mkdir_p, run_, shelly,
                                            silently)
import           System.FilePath
import           System.IO                 (hPutStrLn, stderr)
import           System.IO.Temp            (withTempDirectory)
import           Text.Printf               (printf)
import           Turtle                    hiding (FilePath, format, printf,
                                            stderr)
import qualified Turtle                    as T

import           Bio.Pipeline.Utils        (mapOfFiles)


--------------------------------------------------------------------------------
-- DNA-seq
--------------------------------------------------------------------------------

data BWAOpts = BWAOpts
    { _bwaCores    :: Int          -- ^ number of cpu cores
    , _bwaSeedLen  :: Int     -- ^ seed length, equivalent to -l
    , _bwaMaxMis   :: Int    -- ^ max mismatches in seed, equivalent to -k
    , _bwaReadTrim :: Int       -- ^ dynamic read trimming, equivalent to -q
    , _bwaTmpDir   :: FilePath  -- ^ temp dir
    } deriving (Show)

makeLenses ''BWAOpts

defaultBWAOpts :: BWAOpts
defaultBWAOpts = BWAOpts
    { _bwaCores = 1
    , _bwaSeedLen = 32
    , _bwaMaxMis = 2
    , _bwaReadTrim = 5
    , _bwaTmpDir = "./"
    }

type BWAOptSetter = State BWAOpts ()

{-
-- Determine whether bwa index has been generated
isBWAIndexExist :: FilePath -> IO Bool
isBWAIndexExist dir = do
    fileExist <- testfile (fromText $ T.pack dir ++ "/")
    ls
  where
    filename = "the_suffix_of_this_file_is_the_prefix_of_bwa_indices."
    -}

-- | Generate BWA genome index
bwaMkIndex :: FilePath
           -> FilePath   -- ^ Index prefix, e.g., /path/genome.fa
           -> IO FilePath
bwaMkIndex input prefix = do
    fileExist <- testfile (fromText $ T.pack prefix)
    if fileExist
        then hPutStrLn stderr "BWA index exists. Skipped."
        else do
            mktree $ fromText dir
            cp (fromText $ T.pack input) $ fromText $ T.pack prefix
            hPutStrLn stderr "Generating BWA index"
            shells (T.format ("bwa index -p "%s%" -a bwtsw "%s)
                (T.pack prefix) (T.pack input)) empty
    return prefix
  where
    dir = fst $ T.breakOnEnd "/" $ T.pack prefix

-- | Tag alignment with BWA aligner.
bwaAlign :: (NGS e, IsDNASeq e, Experiment e)
         => FilePath  -- ^ Directory to save the results
         -> FilePath  -- ^ Genome index
         -> BWAOptSetter
         -> e
         -> IO e
bwaAlign dir' index' setter = mapOfFiles fn
  where
    fn e r (Single fl) = if fl^.format == FastqFile || fl^.format == FastqGZip
        then do
            mktree dir
            let output = fromText $ T.format (fp%"/"%s%"_rep"%d%".bam")
                    dir (e^.eid) (r^.number)
                input = fromText $ T.pack $ fl^.location

            stats <- with ( mktempdir (fromText $ T.pack $ opt^.bwaTmpDir) "bwa_align_tmp_dir" ) $ \tmpdir -> do
                let tmp_sai = T.format (fp%"/tmp.sai") tmpdir
                    tmp_sort_bam = T.format (fp%"/sort_bam_tmp") tmpdir
                -- Align reads and save the results to tmp_sai.
                shells ( T.format (
                    "bwa aln -q "%d%" -l "%d%" -k "%d%" -t "%d%" "%fp%" "%fp%" > "%s )
                    (opt^.bwaReadTrim) (opt^.bwaSeedLen) (opt^.bwaMaxMis) (opt^.bwaCores)
                    index input tmp_sai ) empty
                -- Convert sai to sorted bam.
                shells ( T.format (
                    "bwa samse "%fp%" "%s%" "%fp%
                    " | samtools view -Su - | samtools sort - -@ "%d%" -o "%fp%" -T "%s)
                    index tmp_sai input (opt^.bwaCores) output tmp_sort_bam ) empty

            return [ Single $ format .~ BamFile $
                              location .~ T.unpack (T.format fp output) $ fl
                   ]
        else return []
    fn e r (Pair f1 f2) =
        if (f1^.format == FastqFile || f1^.format == FastqGZip) &&
           (f2^.format == FastqFile || f2^.format == FastqGZip)
           then do
                mktree dir
                let output = fromText $ T.format (fp%"/"%s%"_rep"%d%".bam")
                        dir (e^.eid) (r^.number)
                    input1 = fromText $ T.pack $ f1^.location
                    input2 = fromText $ T.pack $ f2^.location

                stats <- with ( mktempdir (fromText $ T.pack $ opt^.bwaTmpDir) "bwa_align_tmp_dir" ) $ \tmpdir -> do
                    let tmp_sort_bam = T.format (fp%"/sort_bam_tmp") tmpdir
                    -- Align reads and save the results to tmp_sai.
                    shells ( T.format (
                        "bwa mem -M -k "%d%" -t "%d%" "%fp%" "%fp%" "%fp%
                        " | samtools view -Su -"%
                        " | samtools sort - -@ "%d%" -o "%fp%" -T "%s )
                        (opt^.bwaSeedLen) (opt^.bwaCores) index input1 input2
                        (opt^.bwaCores) output tmp_sort_bam ) empty

                return [ Single $ format .~ BamFile $
                                  location .~ T.unpack (T.format fp output) $ f1
                       ]
            else return []
    opt = execState setter defaultBWAOpts
    index = fromText $ T.pack index'
    dir = fromText $ T.pack dir'

-- | Remove low quality and redundant tags, fill in mate information.
filterBam :: (NGS e, IsDNASeq e, Experiment e)
          => FilePath  -- ^ directory to save the results
          -> e -> IO e
filterBam dir = mapOfFiles fn
  where
    fn e r (Single fl)
        | fl^.format == BamFile = do
            shelly $ mkdir_p $ fromText $ T.pack dir
            let output = T.pack $ printf "%s/%s_rep%d_filt.bam" dir
                    (T.unpack $ e^.eid) (r^.number)
                input = T.pack $ fl^.location
            bamFilter (pairedEnd e) input output
            return [ Single $ info .~ [] $ format .~ BamFile $
                location .~ T.unpack output $ fl ]
        | otherwise = return []
    fn _ _ _ = return []
    bamFilter isPair input output = withTempDirectory dir "tmp_filt_dir." $ \tmp ->
        shelly $ escaping False $ silently $ do
            let tmp_filt = T.pack $ tmp ++ "/tmp_filt.bam"
                tmp_fixmate = T.pack $ tmp ++ "/tmp_fixmate.bam"
                tmpOutput | isPair = [ "|", "samtools", "sort", "-", "-n"
                    , "-T", T.pack tmp, "-o", tmp_filt ]
                          | otherwise = [">", output]
            run_ "samtools" $ ["view"] ++ outputMode ++
                ["-F", "0x70c", "-q", "30"] ++ compression ++ [input] ++
                tmpOutput
            when isPair $ do
                run_ "samtools" ["fixmate", "-r", tmp_filt, tmp_fixmate]
                run_ "samtools" [ "view", "-F", "1804", "-f", "2", "-u"
                    , tmp_fixmate, "|", "samtools", "sort", "-", "-T", T.pack tmp
                    , "-o", output ]

      where
        outputMode = if isPair then ["-f", "2"] else []
        compression = if isPair then ["-u"] else ["-b"]

-- | Remove duplicates
removeDuplicates :: (NGS e, IsDNASeq e, Experiment e)
                 => FilePath -> FilePath -> e -> IO e
removeDuplicates picardPath dir = mapOfFiles fn
  where
    fn e r (Single fl) = if fl^.format == BamFile
        then do
            shelly $ mkdir_p $ fromText $ T.pack dir
            let output = printf "%s/%s_rep%d_filt_mono.bam" dir
                    (T.unpack $ e^.eid) (r^.number)
                input = fl^.location
                qcFile = printf ("%s/%s_picard.qc") dir (T.unpack $ e^.eid)

            withTempDirectory "./" "tmp_picard_dir." $ \tmp -> shelly $ do
                let markdupTmp = tmp++"/dup_marked.bam"
                    filtTmp = tmp++"/dup_filt.bam"
                -- Mark duplicates
                run_ "java" ["-Xmx4G", "-jar", T.pack picardPath
                    , "MarkDuplicates", T.pack $ "INPUT=" ++ input
                    , T.pack $ "OUTPUT=" ++ markdupTmp
                    , T.pack $ "METRICS_FILE=" ++ qcFile
                    , "VALIDATION_STRINGENCY=LENIENT"
                    , "ASSUME_SORT_ORDER=coordinate", "REMOVE_DUPLICATES=false"]

                -- Remove duplicates.
                escaping False $ run_ "samtools" [ "view", "-F", "0x70c", "-b"
                    , T.pack markdupTmp, ">", T.pack filtTmp ]

                -- Re-sort by names for pairedend sequencing
                if pairedEnd e
                    then run_ "samtools" [ "sort", T.pack filtTmp, "-n"
                        , "-T", T.pack tmp, "-o", T.pack output ]
                    else mv (fromText $ T.pack filtTmp) $ fromText $ T.pack output

                let finalBam = format .~ BamFile $
                               tags .~ ["processed bam file"] $
                               location .~ output $ fl
                    dupQC = format .~ Other $
                            info .~ [] $
                            tags .~ ["picard qc file"] $
                            location .~ qcFile $ fl

                return [Single finalBam, Single dupQC]
        else return []
    fn _ _ _ = return []

bam2Bed :: String -> FileSet -> IO (Maybe FileSet)
bam2Bed prefix (Single fl)
    | fl^.format == BamFile = do
        shelly $ mkdir_p $ fromText $ T.pack $ takeDirectory prefix
        let output = prefix ++ takeBaseName (fl^.location) ++ ".bed.gz"
            bedFile = format .~ BedGZip $
                      location .~ output $ fl
        runBam $ readBam (fl^.location) =$= bamToBed =$= mapC toLine =$=
            unlinesAsciiC =$= gzip $$ sinkFileBS output
        return $ Just $ Single bedFile
    | otherwise = return Nothing
bam2Bed _ _ = return Nothing
{-# INLINE bam2Bed #-}

-- | Convert name sorted BAM to BEDPE suitable for MACS2.
sortedBam2BedPE :: String -> FileSet -> IO (Maybe FileSet)
sortedBam2BedPE prefix (Single fl)
    | fl^.format == BamFile = do
        shelly $ mkdir_p $ fromText $ T.pack $ takeDirectory prefix
        let output = prefix ++ takeBaseName (fl^.location) ++ ".bed.gz"
            bedFile = format .~ BedGZip $
                      location .~ output $ fl
        runBam $ readBam (fl^.location) =$= sortedBamToBedPE =$=
            concatMapC f =$= mapC toLine =$= unlinesAsciiC =$= gzip $$
            sinkFileBS output
        return $ Just $ Single bedFile
    | otherwise = return Nothing
  where
    f (b1, b2)
        | chrom b1 /= chrom b2 || bedStrand b1 == bedStrand b2 = Nothing
        | otherwise =
            let left = if fromJust (bedStrand b1)
                    then chromStart b1 else chromStart b2
                right = if not (fromJust $ bedStrand b2)
                    then chromEnd b2 else chromEnd b1
            in if left <= right
                  then Just $ BED3 (chrom b1) left right
                  else error "Left coordinate is smaller than right coordinate."
sortedBam2BedPE _ _ = return Nothing
{-# INLINE sortedBam2BedPE #-}

-- | Merge multiple gzipped BED files.
mergeReplicatesBed :: Experiment e => FilePath -> e -> IO e
mergeReplicatesBed dir' e = do
    mktree dir
    let fls = map getPath $ e^..replicates.folded.files.folded.filtered isBed
        output = T.format (fp%"/"%s%"_rep0.bed.gz") dir (e^.eid)
        bedFile = Single $ format .~ BedGZip $ location .~ T.unpack output $
            emptyFile
        r = files .~ [bedFile] $ emptyReplicate
    shells (T.format ("zcat "%s%" | gzip -c > "%s)
        (T.unwords $ map T.pack fls) output) empty
    return $ replicates .~ [r] $ e
  where
    dir = fromText $ T.pack dir'
    isBed (Single x) = x^.format == BedGZip
    isBed _ = False
    getPath (Single x) = x^.location


--------------------------------------------------------------------------------
-- RNA-seq
--------------------------------------------------------------------------------

data STAROpts = STAROpts
    { _starCmd    :: T.Text
    , _starCores  :: Int
    , _starTmpDir :: FilePath
    , _starSort   :: Bool
    }

makeLenses ''STAROpts

type STAROptSetter = State STAROpts ()

defaultSTAROpts :: STAROpts
defaultSTAROpts = STAROpts
    { _starCmd = "STAR"
    , _starCores = 1
    , _starTmpDir = "./"
    , _starSort = False
    }


-- | Create index files for STAR
starMkIndex :: FilePath   -- ^ STAR command path
            -> FilePath   -- ^ Directory used to store genome indices
            -> [FilePath] -- ^ Fastq files
            -> FilePath   -- ^ Annotation file
            -> Int        -- ^ The length of the genomic sequence
                          -- around the annotated junction to be used in
                          -- constructing the splice junctions database. Set it
                          -- to "ReadLength-1" or 100 for general purpose.
            -> IO FilePath
starMkIndex star dir fstqs anno r = do
    dirExist <- testdir (fromText $ T.pack dir)
    if dirExist
        then hPutStrLn stderr "STAR index directory exists. Skipped."
        else do
            mktree $ fromText $ T.pack dir
            hPutStrLn stderr "Generating STAR indices"
            shells cmd empty
    return dir
  where
    cmd = T.format (s%" --runThreadN 1 --runMode genomeGenerate --genomeDir "%s%
        " --genomeFastaFiles "%s%" --sjdbGTFfile "%s%" --sjdbOverhang "%d)
        (T.pack star) (T.pack dir) fstqs' (T.pack anno) r
    fstqs' = T.intercalate " " $ map T.pack fstqs

-- | Align RNA-seq raw reads with STAR
starAlign :: FilePath                    -- ^ Output directory
          -> FilePath                    -- ^ STAR genome index
          -> STAROptSetter               -- ^ Options
          -> RNASeq
          -> IO RNASeq
starAlign dir' index' setter = mapOfFiles fn
  where
    isFastq :: File -> Bool
    isFastq x = x^.format == FastqFile || x^.format == FastqGZip
    opt = execState setter defaultSTAROpts
    dir = fromText $ T.pack dir'
    index = fromText $ T.pack index'
    fn e r flset = if null fls then return [] else do
        mktree dir
        with (mktempdir (fromText $ T.pack $ opt^.starTmpDir) "STAR_align_tmp_dir") $ \tmp_dir -> do
            let outputGenome = fromText $ T.format (fp%"/"%s%"_rep"%d%"_genome.bam")
                    dir (e^.eid) (r^.number)
                outputAnno = fromText $ T.format (fp%"/"%s%"_rep"%d%"_anno.bam")
                    dir (e^.eid) (r^.number)
                inputs = T.pack $ unwords $ map (^.location) fls

            shells ( T.format (
                s%" --genomeDir "%fp%" --readFilesIn "%s%
                " --outFileNamePrefix "%fp%"/ --runThreadN "%d%
                (if (head fls)^.format == FastqGZip then " --readFilesCommand zcat" else "")%
                " --genomeLoad NoSharedMemory"%

                " --outFilterType BySJout"%     -- reduces the number of ”spurious” junctions
                " --outFilterMultimapNmax 20"%  -- max number of multiple alignments allowed for a read: if exceeded, the read is considered unmapped
                " --alignSJoverhangMin 8"%      -- minimum overhang for unannotated junctions
                " --alignSJDBoverhangMin 1"%    -- minimum overhang for annotated junctions
                " --outFilterMismatchNmax 999"% -- maximum number of mismatches per pair, large number switches off this filter
                " --outFilterMismatchNoverReadLmax 0.04"% -- max number of mismatches per pair relative to read length: for 2x100b, max number of mismatches is 0.06*200=8 for the paired read
                " --alignIntronMin 20"%         -- minimum intron length
                " --alignIntronMax 1000000"%    -- maximum intron length
                " --alignMatesGapMax 1000000"%  -- maximum genomic distance between mates

                " --outSAMunmapped Within  --outSAMattributes NH HI AS NM MD"%
                " --outSAMheaderCommentFile COfile.txt"%
                " --outSAMheaderHD @HD VN:1.4 SO:coordinate"%
                (if opt^.starSort then " --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 60000000000" else " --outSAMtype BAM Unsorted")%
                (if pairedEnd e then "" else " --outSAMstrandField intronMotif")%
                " --quantMode TranscriptomeSAM --sjdbScore 1"
                ) (opt^.starCmd) index inputs tmp_dir (opt^.starCores) ) empty

            let starOutput | opt^.starSort = "/Aligned.sortedByCoord.out.bam"
                           | otherwise = "/Aligned.out.bam"
            mv (fromText $ T.format (fp%starOutput) tmp_dir) outputGenome

            -- Sorting annotation bam
            if opt^.starSort
                then shells ( T.format ("samtools sort -@ "%d%" -T "%fp%
                        "/sort_bam_tmp -o "%fp%
                        " "%fp%"/Aligned.toTranscriptome.out.bam")
                        (opt^.starCores) tmp_dir outputAnno tmp_dir ) empty
                else mv (fromText $ T.format
                        (fp%"/Aligned.toTranscriptome.out.bam") tmp_dir) outputAnno

            -- Collect bam flagstats
            (_, stats1) <- shellStrict
                (T.format ("samtools flagstat "%fp) outputGenome) empty
            (_, stats2) <- shellStrict
                (T.format ("samtools flagstat "%fp) outputAnno) empty

            let genomeAlignFile = Single $ info .~ [("stat", stats1)] $
                    format .~ BamFile $
                    location .~ T.unpack (T.format fp outputGenome) $
                    tags .~ ["RNA genome align bam"] $ emptyFile
                annoFile = Single $ info .~ [("stat", stats2)] $
                    format .~ BamFile $
                    location .~ T.unpack (T.format fp outputAnno) $
                    tags .~ ["RNA anno align bam"] $ emptyFile
            return [genomeAlignFile, annoFile]
      where
        fls = case flset of
            Single f -> if not (pairedEnd e) && isFastq f then [f] else []
            Pair f1 f2 -> if pairedEnd e && isFastq f1 && isFastq f2
                then [f1,f2] else []


rsemMkIndex :: FilePath   -- ^ Prefix
            -> FilePath   -- ^ annotation file in GFF3 format
            -> [FilePath] -- ^ fastq files
            -> IO FilePath
rsemMkIndex prefix anno fstqs = do
    dirExist <- testdir (fromText $ T.pack dir)
    if dirExist
        then hPutStrLn stderr "RSEM index directory exists. Skipped."
        else do
            mktree $ fromText $ T.pack dir
            hPutStrLn stderr "Generating RSEM indices"
            shells cmd empty
    return prefix
  where
    cmd = T.format ("rsem-prepare-reference --gtf "%s%" "%s%" "%s)
        (T.pack anno) (T.intercalate "," $ map T.pack fstqs) (T.pack prefix)
    dir = T.unpack $ fst $ T.breakOnEnd "/" $ T.pack prefix


data RSEMOpts = RSEMOpts
    { _rsemPath  :: T.Text
    , _rsemCores :: Int
    , _rsemSeed  :: Int32
    }

makeLenses ''RSEMOpts

type RSEMOptSetter = State RSEMOpts ()

defaultRSEMOpts :: RSEMOpts
defaultRSEMOpts = RSEMOpts
    { _rsemPath = ""
    , _rsemCores = 1
    , _rsemSeed = 12345
    }

-- | Gene and transcript quantification using rsem
rsemQuant :: FilePath
          -> FilePath
          -> RSEMOptSetter
          -> RNASeq
          -> IO RNASeq
rsemQuant dir' indexPrefix setter = mapOfFiles fn
  where
    opt = execState setter defaultRSEMOpts
    dir = fromText $ T.pack dir'
    isAnnoBam f = f^.format == BamFile && "RNA anno align bam" `elem` f^.tags
    fn e r (Single fl) = if not (isAnnoBam fl) then return [] else do
        mktree dir
        let input = fromText $ T.pack $ fl^.location
            output = T.format (fp%"/"%s%"_rep"%d%"_rsem") dir (e^.eid) (r^.number)
        shells ( T.format ( s%
            "rsem-calculate-expression --bam --estimate-rspd --calc-ci"%
            " --seed "%d%" -p "%d%" --no-bam-output --ci-memory 30000"%
            (if pairedEnd e then " --paired-end --forward-prob 0 " else " ")%fp%
            " "%s%" "%s
            ) (opt^.rsemPath) (opt^.rsemSeed) (opt^.rsemCores) input
            (T.pack indexPrefix) output ) empty

        let geneQuant = Single $ format .~ Other $
                location .~ T.unpack output ++ ".genes.results" $
                tags .~ ["gene quantification"] $ fl
            transcirptQuant = Single $ format .~ Other $
                location .~ T.unpack output ++ ".isoforms.results" $
                tags .~ ["transcript quantification"] $ fl
        return [geneQuant, transcirptQuant]
    fn _ _ _ = return []
