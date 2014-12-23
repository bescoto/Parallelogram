\name{Parallelogram}
\alias{Parallelogram}
\title{Create a Parallelogram object}

\description{ Return a Parallelogram object, which can be used to
  compute raw/on-level written/earned premium by time period. }

\usage{
Parallelogram(rate.change.df, term.len, written.df=NULL)
}

\arguments{
  \item{rate.change.df}{a data frame containing historical rate change
  information.  It should have columns \code{year} and
  \code{rate.change} containing the dates and magnitudes of historical
  rate changes.}

  \item{term.len}{The length of a policy term in years}

  \item{written.df}{Optional data frame about rate of written premium.
    It should contain the columns \code{year} and \code{written},
    meaning that starting on \code{year}, the rate of written premium
    was \code{written}.  If no argument is given, assume constant
    written premium rate of 1.
}}

\value{
  A parallelogram object.
  }

\details{
  This implements the classic "parallelogram method" that is familiar to
  P&C actuaries.  However, note that:

  \itemize{
  \item Time is modeled with floating point numbers.  If \code{Date}
  objects are given, they will be turned into a float.  This means that
  leap years are not distinguished from ordinary years, and, for
  instance, a term length of .5 will often not correspond to the next 6
  calendar months.

  \item The rate at which premium is written is allowed to vary as a
  step function.  If written premium is fixed the classic
  "parallelogram" shape is obtained.
  
}
  
  This doesn't check for exotic data frame column types (e.g. one column
  with two sub-columns), and may not handle them correctly).
}

\seealso{
  \code{\link{EarnedByPeriod}},
  \code{\link{WrittenByPeriod}}
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

# Check a few values
written.df <- WrittenByPeriod(p, periods.out)
stopifnot(all(written.df$raw.premium == 1))
stopifnot(with(written.df, all(onlevel.premium == rate.level)))
stopifnot(all.equal(written.df$onlevel.premium,
                    c(1.304266656, 1.2189408, 1.2377904, 1.12226399999999, 1.05)))

earned.df <- EarnedByPeriod(p, periods.out)
stopifnot(all(earned.df$raw.premium == 1))
stopifnot(with(earned.df, all(onlevel.premium == rate.level)))
stopifnot(all.equal(earned.df$onlevel.premium,
                    c(1.304266656, 1.261603728, 1.2236532, 1.19343839999999, 1.0899332)))


##### Example 2: Same rate changes as Ex1, but variable rate of written premium

# variable written premium
written.df <- data.frame(year=2001:2005,
                         written=c(30, 45, 75, 30, 55))
p <- Parallelogram(rate.change.df, term.len, written.df)
plot(p)

EarnedByPeriod(p, periods.out)
WrittenByPeriod(p, periods.out)


# Check a few values
written.df <- WrittenByPeriod(p, periods.out)
stopifnot(all(written.df$raw.premium == c(30, 45, 75, 30, 55)))
stopifnot(with(written.df, all(onlevel.premium == raw.premium * rate.level)))
stopifnot(all.equal(written.df$rate.level,
                    c(1.304266656, 1.2189408, 1.2377904, 1.12226399999999, 1.05)))

earned.df <- EarnedByPeriod(p, periods.out)
stopifnot(all.equal(earned.df$raw.premium, c(30, 37.5, 60, 52.5, 42.5)))
stopifnot(with(earned.df, all.equal(onlevel.premium, raw.premium * rate.level)))
stopifnot(all.equal(earned.df$rate.level,
                    c(1.304266656, 1.2530711424, 1.2248313, 1.21648559999999, 1.08554108235294)))

}
