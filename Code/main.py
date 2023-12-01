# python3 ./Code/main.py 
try:
    from pip import main as pipmain
except ImportError:
    from pip._internal import main as pipmain

from performance_tester import performance_test

if __name__ == "__main__":
    pipmain(['install', '-r', 'requirements.txt'])
    performance_test()