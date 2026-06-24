docs:
	julia --project=docs docs/make.jl

check-math: docs
	node docs/check-math.js

normalize-quotes:
	python3 docs/normalize-quotes.py

.PHONY: docs check-math normalize-quotes
