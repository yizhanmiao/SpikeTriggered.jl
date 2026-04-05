.PHONY: all, docs, update, serve

all: update docs

update:
	julia --project=@. -e 'using Pkg; Pkg.update()'

docs: clean
	julia --project=@. make.jl

clean:
	-rm src/CHANGELOG.md
	-rm -r src/examples/*.md
	-rm -r build

serve:
	python3 -m http.server
