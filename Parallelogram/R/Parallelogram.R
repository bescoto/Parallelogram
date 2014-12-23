#######################################################################
# All code for the Parallelogram package


Parallelogram <- function(rate.change.df, term.len, written.df=NULL) {
  # Return on-level written and earned premium by interval
  #
  # Args:
  # rate.change.df - a data frame with dates and rate changes.  Fields
  #   are year and rate.change
  # term.len - length of policy term in years
  # written.df - a data frame with written premium by period.  Fields are
  #   year and written.
  written.df <- .CheckWrittenDF(written.df, term.len)
  rate.level.df <- .RateLevelDF(rate.change.df)
  written.steps <- .WrittenSteps(written.df, rate.level.df)
  inforce.pieces <- .InforcePieces(written.steps, term.len)
  inforce.funcs <- lapply(inforce.pieces, function(pair) pair$f)
  written.knots <- lapply(written.steps, function(step) knots(step))
  inforce.knots <- lapply(inforce.pieces, function(pair) pair$knots)
  result <- list(term.len=term.len,
                 rate.level.df=rate.level.df,
                 written.df=written.df,
                 written.funcs=written.steps,
                 written.knots=written.knots,
                 written.knots.all=Reduce(function(x, y) c(x, y),
                     written.knots, NULL),
                 TotalWritten=.SumFunctions(written.steps),
                 inforce.funcs=inforce.funcs,
                 inforce.knots=inforce.knots,
                 inforce.knots.all=Reduce(function(x, y) c(x, y),
                     inforce.knots, NULL),
                 TotalInforce=.SumFunctions(inforce.funcs))
  class(result) <- "parallelogram"
  return(result)
}

.Assert <- function(boolean, error.msg) {
  # Abort with error message if boolean not true
  if (!boolean)
    stop(error.msg)
}

.CheckWrittenDF <- function(written.df, term.len) {
  # Make sure written.df is in the right format
  if (is.null(written.df)) # Default to constant writing
    return(data.frame(written=1, year=0))
  .Assert("data.frame" %in% class(written.df),
         "written.df should be a data frame with premium in it")
  .Assert(nrow(written.df) >= 1, "written.df has no rows!")
  written.df$year <- .CheckYearField(written.df$year)
  .Assert("written" %in% names(written.df) || "earned" %in% names(written.df),
          'written.df data frame requires either "written" or "earned" column')
  if (!("written" %in% names(written.df)))
    written.df$written <- written.df$earned - term.len / 2
  .Assert(!any(is.na(written.df$written)) && !any(is.na(written.df$year)),
          "NAs not allowed inside written.df")
  n <- nrow(written.df)

  if (n > 1)
    .Assert(all(written.df$year[2:n] > written.df$year[1:(n - 1)]),
            "years in written.df must be in order")
  return(written.df)
}

.CheckYearField <- function(year) {
  # Check the type of the year field, converting if necessary
  if ("Date" %in% class(year)) {
    warning("Date field found, converting to decimal year")
    return(lubridate::decimal_date(year))
  }
  .Assert(any(c("numeric", "integer") %in% class(year)),
          "Year must be specified as decimal, e.g. 2012.34")
  return(year)
}

.RateLevelDF <- function(rate.change.df) {
  # Check rate.change.df and return rate.level.df
  #
  # rate.level.df will have one row for every different rate level.
  # The columns are:
  #    start - starting time of the period, possibly -Inf
  #    end - end time of the period, possibly Inf
  #    rate.level - the rate level, where 1.0 is the current rate level
  if (is.null(rate.change.df) || nrow(rate.change.df) == 0)
    return(data.frame(start=-Inf, end=Inf, rate.level=1))
  .Assert(all(c("year", "rate.change") %in% names(rate.change.df)),
         'rate.change.df must contain "year" and "rate.change" columns')
  rate.change.df$year <- .CheckYearField(rate.change.df$year)
  rate.change.df <- rate.change.df[order(rate.change.df$year), ]
  start <- c(-Inf, rate.change.df$year)
  end <- c(rate.change.df$year, Inf)
  rate.level <- c(rev(cumprod(rev(1 + rate.change.df$rate.change))), 1.0)
  return(data.frame(start=start, end=end, rate.level=rate.level))
}

.WrittenSteps <- function(written.df, rate.level.df) {
  # Return list of step functions, one for each rate change period
  # plus additional function for total written
  Helper <- function(start, end) {
    # Return a single step function covering written from start to end
    prev.written <- if (any(written.df$year < start))
      written.df[max(which(written.df$year < start)), "written"]
      else written.df$written[1]
    sub.df <- written.df[start <= written.df$year & written.df$year < end, ]
    return(stepfun(x=c(start, sub.df$year, end),
                   y=c(0, prev.written, sub.df$written, 0)))

  }
  return(mlply(rate.level.df[, c("start", "end")], Helper))
}

.InforcePieces <- function(written.steps, term.len) {
  # Return list of inforce information, one element for each rating period
  #
  # Each element in the result be a list with two subelements:
  #   f - a piecewise linear functions giving inforce premium at that time
  #   knots - a numeric vector of f's turning points
  Helper <- function(stepfun) {
    # Return f and knots for a single step function
    f.knots <- GetKnots(stepfun)
    if (length(f.knots) == 0) # written and inforce are constant
      return(list(f=function(x) stepfun(0) * term.len, knots=NULL))
    y <- YVals(f.knots, stepfun)
    f <- approxfun(f.knots, y, rule=2, method="linear")
    return(list(knots=f.knots, f=f))
  }

  GetKnots <- function(stepfun) {
    # Return the knots of the inforce premium function given written step fun
    f.knots <- sort(unique(c(knots(stepfun), knots(stepfun) + term.len)))
    f.knots <- f.knots[-Inf < f.knots & f.knots < Inf]
    .Assert(length(f.knots) >= 2, "Sanity check--this shouldn't be false")
    return(f.knots)
  }

  YVals <- function(f.knots, stepfun) {
    # Given knots of inforce function and stepfun, return inforce prem function
    #
    # This integrates the stepfun by rolling over the interval (x -
    # term.len, x) one unit at a time (once per iteration of the for
    # loop).
    x <- f.knots[1]
    y <- stepfun(x - 1) * term.len # stepfun is constant before first knot
    for (new.x in f.knots[2:length(f.knots)]) {
      new.inforce <- stepfun(tail(x, 1)) * (new.x - tail(x, 1))
      expired <- stepfun(tail(x, 1) - term.len) * (new.x - tail(x, 1))
      x <- c(x, new.x)
      y <- c(y, tail(y, 1) + new.inforce - expired)
    }
    return(y)
  }

  return(lapply(written.steps, Helper))
}

.SumFunctions <- function(func.list) {
  # Return a function yielding the sum of the given functions
  ReturnFunc <- function(vals) {
    # Return total written for each year specified
    result <- rep(0, length(vals))
    for (f in func.list)
      result <- result + f(vals)
    return(result)
  }
}

plot.parallelogram <- function(x, xrange=NULL,
                               fill.colors=c("#7C6C5C", "#583008"),
                               text.color="#06401B",
                               written.prem.color="#583008",
                               ...) {
  # Plot the different rating periods embedded in the parallelogram
  xrange <- .GetXRange(x, xrange)
  graph.df <- .MakeInforceDF(x, xrange)
  period.colors <- rep(fill.colors, nrow(x$rate.level.df))
  names(period.colors) <- NULL # otherwise the names confuse ggplot

  plot <- (ggplot(data=graph.df)
           + geom_area(aes(x=year, y=inforce, group=rating.period,
                         fill=rating.period),
                     alpha=0.5, position="stack")
           + scale_fill_manual(name="Rating Period", values=period.colors,
                               guide=FALSE)
           + labs(x="Year", y="Inforce Premium"))

  plot <- .PlotAddText(x, plot, graph.df, text.color)
  plot <- .PlotAddWritten(x, plot, graph.df, written.prem.color)
  return(plot)
}

.PlotAddText <- function(p, plot, graph.df, text.color) {
  # Return plot including ggplot text
  midpoint.df <- .FindPeriodMidpoints(p, graph.df)
  midpoint.df$rate.level <- formatC(p$rate.level.df$rate.level,
                                    digits=3, format="f")
  return(plot + geom_text(data=midpoint.df, color=text.color, size=4,
                          aes(x=x, y=y, label=rate.level)))
}

.PlotAddWritten <- function(p, plot, graph.df, line.color) {
  # add a line for written premium
  graph.df$written <- p$TotalWritten(graph.df$year)
  plot <- (plot + geom_line(aes(x=year, y=written),
                            data=graph.df, color=line.color, linetype=2))
  return(plot)
}

.GetXRange <- function(p, xrange) {
  # Return a valid xrange for the graph
  if (!is.null(xrange))
    return(xrange)

  if (nrow(p$rate.level.df) > 1)
    return(with(p$rate.level.df,
                c(min(end) - .1 * p$term.len, max(start) + 1.1 * p$term.len)))
}

.MakeInforceDF <- function(p, xrange, spacing=0.02) {
  # Return a data frame of inforce premium by rating period, suitable for ggplot
  xvals <- seq(from=xrange[1], to=xrange[2], by=spacing)
  inforce.graph.df <- data.frame(rating.period=NULL, year=NULL, inforce=NULL)
  for (i in seq(along=p$inforce.funcs)) { # add points for rating period i - 1
    new.df <- data.frame(rating.period=i, year=xvals,
                         inforce=p$inforce.funcs[[i]](xvals))
    inforce.graph.df <- rbind(inforce.graph.df, new.df)
  }
  inforce.graph.df$rating.period <- as.factor(inforce.graph.df$rating.period)
  # now reverse so parallelograms plot in the right direction
  return(inforce.graph.df[nrow(inforce.graph.df):1, ])
}

.FindPeriodMidpoints <- function(p, inforce.graph.df) {
  # Return the midpoints of the rating period, used for text placement
  Helper <- function(df) {
    # Return midpoint of one rating period
    x <- with(df, sum(year * inforce) / sum(inforce))
    y <- p$TotalInforce(x) / 2
    return(data.frame(x=x, y=y))
  }
  result.df <- ddply(inforce.graph.df, .(rating.period), Helper)

  # Now jitter the positions to prevent text overlap
  jitter <- ifelse(seq(length.out=nrow(result.df)) %% 2 == 1, 1, -1)
  result.df$x <- result.df$x + 0.08 * term.len * jitter
  result.df$y <- result.df$y * (1 + 2 * 0.08 * jitter)
  return(result.df)
}

.EarnedByPeriodDetails <- function(parallelogram, periods.out) {
  # Return earned premium by period from parallelogram results
  #
  # Inputs:
  #   parallelogram - parallelogram results object
  #   periods.out - vector of period begin/end numbers
  # Output will be a data frame with these columns:
  #   start, end - the year of the start and end of the period.out
  #   rating.period - the number of the rating period
  #   earned - the premium earned in that out period that was written under
  #     that rating period.
  Helper <- function(period.df) {
    # Given start and end times for the period, return section of final result
    .Assert(nrow(period.df) == 1, "Two periods with same start and end")
    RatingHelper <- function(rating.period) {
      # Returned earned premium during the period for a single rating period
      #
      # To do this we integrate over the piecewise linear inforce
      # function.  Since we know where the knots are, to find the area
      # we can just average the beginning and end of each section and
      # multiply by the width (each section is trapezoidal).
      knots <- parallelogram$inforce.knots[[rating.period]]
      knots <- knots[period.df$year.start < knots & knots < period.df$year.end]
      Inforce <- parallelogram$inforce.funcs[[rating.period]]
      oldx <- period.df$year.start
      area <- 0
      for (newx in c(knots, period.df$year.end)) {
        area <- area + (Inforce(oldx) + Inforce(newx)) / 2 * (newx - oldx)
        oldx <- newx
      }
      return(area / parallelogram$term.len)
    }

    rating.periods <- seq(length.out=nrow(parallelogram$rate.level.df))
    earned <- sapply(rating.periods, RatingHelper)
    return(data.frame(rating.period=rating.periods, earned=earned))
  }

  .Assert("parallelogram" %in% class(parallelogram),
          "First parameter should be a parallelogram as made by Parallelogram")
  periods.df <- .CheckPeriodsOut(periods.out)
  return(ddply(periods.df, .(year.start, year.end), Helper))
}

.CheckPeriodsOut <- function(periods.out) {
  # Make sure periods.out is in right format, return periods.df
  n <- length(periods.out)
  .Assert(n >= 2, "periods.out requires at least beginning and end points")
  .Assert(all(periods.out[2:n] > periods.out[1:(n - 1)]),
         "times in periods.out need to be in order")
  .Assert(!any(is.na(periods.out)), "NA's not allowed in periods.out")
  return(data.frame(period.out=1:(n - 1),
                    year.start=periods.out[1:(n - 1)],
                    year.end=periods.out[2:n]))
}

EarnedByPeriod <- function(p, periods.out) {
  # Return data frame summarizing earned premium by period
  .Assert("parallelogram" %in% class(p))
  details.df <- .EarnedByPeriodDetails(p, periods.out)
  details.df$onlevel.prem <- with(details.df,
               earned * p$rate.level.df$rate.level[details.df$rating.period])
  result.df <- ddply(details.df, .(year.start, year.end), summarize,
                     raw.premium=sum(earned), onlevel.premium=sum(onlevel.prem))
  result.df$rate.level <- with(result.df, onlevel.premium / raw.premium)
  return(result.df)
}

WrittenByPeriod <- function(p, periods.out) {
  # Return data frame summarizing written premium by period
  .Assert("parallelogram" %in% class(p),
          "First argument to WrittenByPeriod should be a parallelogram object")
  ByInputPeriod <- function(period.df) {
    # Process a single period of periods.out
    .Assert(nrow(period.df) == 1, "Two periods with same start and end")
    df <- transform(p$rate.level.df,
                    new.start=pmax(period.df$year.start, start),
                    new.end=pmin(period.df$year.end, end))
    df <- transform(df,
                    length=pmax(0, new.end - new.start),
                    written.rate=p$TotalWritten(new.start))
    return(with(df, data.frame(
        raw.premium=sum(written.rate * length),
        onlevel.premium=sum(written.rate * length * rate.level))))
  }
  periods.df <- .CheckPeriodsOut(periods.out)
  result.df <- ddply(periods.df, .(year.start, year.end), ByInputPeriod)
  result.df$rate.level = with(result.df, onlevel.premium / raw.premium)
  return(result.df)
}
