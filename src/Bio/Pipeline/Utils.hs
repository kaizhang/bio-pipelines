{-# LANGUAGE OverloadedStrings #-}

module Bio.Pipeline.Utils where

import qualified Data.Text as T

-- | Get the prefix of a filename, e.g., "path/xxx.txt" -> "path/xxx"
getPrefix :: T.Text -> T.Text
getPrefix x = if suffix == "gz"
    then T.init $ fst $ T.breakOnEnd "." $ T.init prefix
    else T.init prefix
  where
    (prefix, suffix) = T.breakOnEnd "." x
