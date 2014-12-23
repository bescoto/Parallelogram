\name{WrittenByPeriod}
\alias{EarnedByPeriod}
\alias{WrittenByPeriod}
\title{Compute Earned or Written Premium by Period}

\description{ Given a parallelogram object, return a data frame
  containing rate level and onlevel premium information by period }

\usage{
EarnedByPeriod(p, periods.out)
WrittenByPeriod(p, periods.out)
}

\arguments{
  \item{p}{A parallelogram object, as produced by the
  \code{\link{Parallelogram}} function}

  \item{periods.out}{A vector of years indicating which periods should
  be computed}
}

\value{
  A data frame with these columns:
  \itemize{
  \item \code{year.start} - the beginning of each time period, as specified
  by \code{periods.out}
  \item \code{year.end} - the end of each time period, as specified by
  \code{periods.out}
  \item \code{raw.premium} - the amount of premium written or earned in
  the period
  \item \code{onlevel.premium} - the onlevel premium written or earned
  in the period
  \item \code{rate.level} - \code{onlevel.premium} divided by \code{raw.premium}
 }

  \code{WrittenByPeriod} returns a data frame about written premium.
  \code{EarnedByPeriod} returns a data frame about earned premium.
}

\details{
  \code{EarnedByPeriod} implements the parallelogram method.  See
  \code{\link{Parallelogram}} for more information.
}

\seealso{
  \code{\link{Parallelogram}}
}

\examples{


##### Example 1: constant rate of written premium (classic parallelogram shape)

rate.change.df <- data.frame(year=c(2002, 2003.5, 2004.1, 2004.4, 2005.5),
                             rate.change=c(.07, -.03, .12, .02, .1))
term.len <- 1
p <- Parallelogram(rate.change.df, term.len)
plot(p)

periods.out <- 2001:2006
WrittenByPeriod(p, periods.out)
EarnedByPeriod(p, periods.out)


##### Example 2: Same rate changes as Ex1, but variable rate of written premium

# variable written premium
written.df <- data.frame(year=2001:2005,
                         written=c(30, 45, 75, 30, 55))
p <- Parallelogram(rate.change.df, term.len, written.df)
plot(p)

EarnedByPeriod(p, periods.out)
WrittenByPeriod(p, periods.out)

}
