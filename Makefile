.PHONY: coverage-report clean clean-all test

test: clean
	pixi run julia --project=@. -e 'using Pkg; Pkg.test(coverage=false)'

coverage-report: clean-all
	pixi run julia --project=@. --code-coverage=user -e 'using Pkg; Pkg.test(coverage=true)'
	julia -e 'using Coverage;cov=process_folder();LCOV.writefile("coverage.info", cov)'
	# julia -e 'using Coverage;cov=process_folder();append!(cov, process_folder("ext"));LCOV.writefile("coverage.info", cov)'
	genhtml coverage.info -o coverage-report
	open coverage-report/index.html

clean:
	-find src -name "*.cov" -type f -delete
	-find test -name "*.cov" -type f -delete
	# -find ext -name "*.cov" -type f -delete # remove coverage files from `ext`
	-find test -type d -name "jl_*" -exec rm -rf {} +
	-rm coverage.info

clean-coverage-report:
	-rm -r coverage-report/

clean-all: clean clean-coverage-report
