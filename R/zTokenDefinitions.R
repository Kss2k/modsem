LavToken        <- setClass("LavToken")
setMethod("assignSubClass", "LavToken", assignSubClass.LavToken)
setMethod("evalToken", "LavToken", evalToken.LavToken)

  LavName         <- setClass("LavName",
                              contains = "LavToken")
  setMethod("assignSubClass", "LavName", assignSubClass.LavName)
  LavObject       <- setClass("LavObject", contains = "LavName")
  LavFunciton     <- setClass("LavFunction", contains = "LavName")
  setMethod("evalToken", "LavFunction", evalToken.LavFunction)

  LavOperator     <- setClass("LavOperator", contains = "LavToken")
  setMethod("assignSubClass", "LavOperator", assignSubClass.LavOperator)
  setMethod("evalToken", "LavOperator", evalToken.LavOperator)


    LavMeasure      <- setClass("LavMeasure", contains = "LavOperator" )
    LavPredict      <- setClass("LavPredict", contains = "LavOperator")
    LavCovar        <- setClass("LavCovar" , contains = "LavOperator")
    LavAdd          <- setClass("LavAdd"   , contains = "LavOperator")
    setMethod("evalToken", "LavAdd", evalToken.LavAdd)
    LavModify       <- setClass("LavModify")
    setMethod("evalToken", "LavModify", evalToken.LavModify)
    LavLessLeft     <- setClass("LavLessLeft" , contains = "LavOperator")
    LavLessRight    <- setClass("LavLessRight", contains = "LavOperator")
    LavEquals       <- setClass("LavEquals"   , contains = "LavOperator")
    LavInteraction  <- setClass("LavInteraction", contains = "LavOperator")
    setMethod("evalToken", "LavInteraction", evalToken.LavInteraction)
    LavSeperator    <- setClass("LavSeperator")
    setMethod("evalToken", "LavSeperator", evalToken.LavSeperator)

  LavClosure      <- setClass("LavClosure", contains = "LavToken")
  setMethod("assignSubClass", "LavClosure", assignSubClass.LavClosure)
    LeftBracket     <- setClass("LeftBracket", contains = "LavClosure")
    setMethod("evalToken", "LeftBracket", evalToken.LeftBracket)
    RightBracket    <- setClass("RightBracket", contains = "LavClosure")
    setMethod("evalToken", "RightBracket", evalToken.RightBracket)

  LavComment      <- setClass("LavComment", contains = "LavToken")
  setMethod("evalToken", "LavComment", evalToken.LavComment)

  LavBlank        <- setClass("LavBlank", contains = "LavToken")
  setMethod("evalToken", "LavBlank", evalToken.LavBlank)

  LavNumeric      <- setClass("LavNumeric", contains = "LavToken")
  setMethod("assignSubClass", "LavNumeric", assignSubClass.LavNumeric)
  LavString       <- setClass("LavString", contains = "LavToken")

  LavMathExpr       <- setClass("LavMathExpr", contains = "LavToken")
  setMethod("assignSubClass", "LavMathExpr", assignSubClass.LavMathExpr)
