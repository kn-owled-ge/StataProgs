capture program drop cls
program define cls, nclass
	version 10.1
	qui query
	qui loc lines = c(pagesize)
	if c(more) == "on" {
		qui set more off
		forvalues i=1/100 {
			display _newline(`lines')
		}
		qui set more on
	}
	else {
		forvalues i=1/100 {
			display _newline(`lines')
		}
	}
end
