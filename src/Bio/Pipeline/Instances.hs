{-# LANGUAGE DeriveGeneric      #-}
{-# LANGUAGE FlexibleInstances  #-}
{-# LANGUAGE StandaloneDeriving #-}
{-# LANGUAGE TemplateHaskell    #-}

module Bio.Pipeline.Instances where

import           Bio.Data.Bed
import           Bio.Motif             (Motif (..), PWM (..))
import           Control.Monad
import           Data.Aeson
import           Data.Aeson.TH
import qualified Data.ByteString.Char8 as B
import           Data.Hashable
import           Data.Maybe
import           Data.Serialize        (Serialize)
import           Data.Binary (Binary)
import qualified Data.Text             as T
import           GHC.Generics          (Generic)

instance FromJSON B.ByteString where
    parseJSON (String x) = return $ B.pack $ T.unpack x
    parseJSON _ = mzero

instance ToJSON B.ByteString where
    toJSON x = String $ T.pack $ B.unpack x

$(deriveJSON defaultOptions ''BED)
$(deriveJSON defaultOptions ''BED3)
$(deriveJSON defaultOptions ''NarrowPeak)

deriving instance Generic BED
instance Serialize BED
instance Binary BED
