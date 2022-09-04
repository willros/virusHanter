install:
	pip install --upgrade pip
	pip install -r requirements.txt
    
    
format:
	black *.py plotting/*.py ldap-login/*.py forms/*.py helper_functions/*.py

lint:
	pylint --disable=R,C *.py plotting/*.py ldap-login/*.py forms/*.py

#test:
#	python -m pytest -vv --cov=virusHanter --cov=main tests/*
#build:
#	docker build -t virushanter .
#run_docker:
#	docker run -p 127.0.0.1:8080:8080 virushanter
#deploy:
#	#deploy
#all: install format lint test build deploy