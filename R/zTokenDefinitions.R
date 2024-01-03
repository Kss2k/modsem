LavToken        <- setClass("LavToken")
setMethod("assignSubClass", "LavToken", assignSubClass.LavToken)
setMethod("evalToken", "LavToken", evalToken.LavToken)

  LavName         <- setClass("LavName")
  setMethod("assignSubClass", "LavName", assignSubClass.LavName)
  LavObject       <- setClass("LavObject")
  LavFunciton     <- setClass("LavFunction")
  setMethod("evalToken", "LavFunction", evalToken.LavFunction)

  LavOperator     <- setClass("LavOperator")
  setMethod("assignSubClass", "LavOperator", assignSubClass.LavOperator)
  setMethod("evalToken", "LavOperator", evalToken.LavOperator)


    LavMeasure      <- setClass("LavMeasure")
    LavPredict      <- setClass("LavPredict")
    LavCovar        <- setClass("LavCovar")
    LavAdd          <- setClass("LavAdd")
    setMethod("evalToken", "LavAdd", evalToken.LavAdd)
    LavModify       <- setClass("LavModify")
    setMethod("evalToken", "LavModify", evalToken.LavModify)
    LavLessLeft     <- setClass("LavLessLeft")
    LavLessRight    <- setClass("LavLessRight")
    LavEquals       <- setClass("LavEquals")
    LavInteraction  <- setClass("LavInteraction")
    setMethod("evalToken", "LavInteraction", evalToken.LavInteraction)
    LavSeperator    <- setClass("LavSeperator")
    setMethod("evalToken", "LavSeperator", evalToken.LavSeperator)

  LavClosure      <- setClass("LavClosure")
  setMethod("assignSubClass", "LavClosure", assignSubClass.LavClosure)
    LeftBracket     <- setClass("LeftBracket")
    setMethod("evalToken", "LeftBracket", evalToken.LeftBracket)
    RightBracket    <- setClass("RightBracket")
    setMethod("evalToken", "RightBracket", evalToken.RightBracket)

  LavComment      <- setClass("LavComment")
  setMethod("evalToken", "LavComment", evalToken.LavComment)

  LavBlank        <- setClass("LavBlank")
  setMethod("evalToken", "LavBlank", evalToken.LavBlank)

  LavNumeric      <- setClass("LavNumeric")

  LavString       <- setClass("LavString")
