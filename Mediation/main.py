from notebooks.1CleanAndCombine import clean_and_combine
from notebooks.1CleanAndCombine import run_clean_and_combine

def main():
    data_cleaner = clean_and_combine(threshold=5, max_missing=384370)
    data_cleaner.run()

    run_clean = run_clean_and_combine(threshold=5, max_missing=384370)
    run_clean.run_clean_and_combine()

if __name__ == "__main__":
    main()