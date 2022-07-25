// Borrrowed from _iv_parse.ado, with modification
/* The syntax for rdhonest is: 
fuzzy RD: Y (X=Z) W
RD: Y X W
where
Y: depvar, X: runvar, Z: treat, W: covar
*/

program rdparse, sclass

	//syntax cmd_to_parse [if] [in]
	//noi di "`cmd_to_parse'"
	//marksample touse
	local n 0

	gettoken lhs 0 : 0, parse(" ,[") match(paren) bind
	if (strpos("(",`"`lhs'"')) {
		fvunab lhs : `lhs'
		if `:list sizeof lhs' > 1 {
			gettoken lhs rest : lhs
			local 0 `"`rest' `0'"'
		}
	}

	IsStop `lhs'
	if `s(stop)' {
		error 198 
	}

	_fv_check_depvar `lhs'
	while `s(stop)'==0 {
		if "`paren'"=="(" {
			local n = `n' + 1
			if `n'>1 {
				capture noi error 198
				di as error `"syntax is "(one running variable = one treatment variable)""'
				exit 198
			}
			gettoken p lhs : lhs, parse(" =") bind
			while "`p'"!="=" {
				if "`p'"=="" {
					capture noi error 198
					di as error `"syntax is "(one running variable = one treatment variable)""'
					di as error `"the equal sign "=" is required"'
					exit 198
				}
				local runvar`n' `runvar`n'' `p'
				gettoken p lhs : lhs, parse(" =") bind
			}
			if "`runvar`n''" != "" {
				fvunab runvar`n' : `runvar`n''
				// * unnote the following lines if we want to limit the number of running variables to 1
				//if `:list sizeof runvar`n'' > 1 {
				//capture noi error 198
				//di as error `"syntax is "(one running variable = one treatment variable)""'
				//di as error `"only one running variable allowed"'
				//}
			}
			fvunab covar`n' : `lhs'
			// * unnote the following lines if we want to limit the number of treatment variables to 1
			//if `:list sizeof covar`n'' > 1 {
			//capture noi error 198
			//di as error `"syntax is "(one running variable = one treatment variable)""'
			//di as error `"only one treatment variable allowed"'
			//}
		}
		else {
			local covar `covar' `lhs'
		}
		gettoken lhs 0 : 0, parse(" ,[") match(paren) bind
		IsStop `lhs'
	}
	mata: st_local("0",strtrim(st_local("lhs")+ " " + st_local("0")))

	fvunab covar : `covar'
	fvexpand `covar'
	local covar `r(varlist)'
	tokenize `covar'
	local lhs "`1'"
	local 1 " "
	local covar `*'
	
	// Eliminate vars from `covar1' that are in `covar'
	local treat : list covar1 - covar
	if ("`runvar1'" != "") {
		fvunab runvar1 : `runvar1'
		fvexpand `runvar1'
		local runvar `r(varlist)'
	}
	else {
		tokenize `covar'
		local runvar "`1'"
		local 1 " "
		local covar `*'
	}
	
	// `lhs' contains the dependent variable
	// `runvar' contains the running variable
	// `treat' contains the treatment variable
	// `covar' contains RHS covariates
	// `0' contains whatever is left over (if/in, weights, options)
	
	sret local depvar `lhs'
	sret local covar `covar'
	sret local runvar `runvar'
	sret local treat `treat'
	sret local zero `"`0'"'

end

// Borrowed from ivregress.ado	
program define IsStop, sclass

	if `"`0'"' == "[" {
		sret local stop 1
		exit
	}
	if `"`0'"' == "," {
		sret local stop 1
		exit
	}
	if `"`0'"' == "if" {
		sret local stop 1
		exit
	}
	if `"`0'"' == "in" {
		sret local stop 1
		exit
	}
	if `"`0'"' == "" {
		sret local stop 1
		exit
	}
	else {
		sret local stop 0
	}

end
