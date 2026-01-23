import argparse
import re


def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="clean HTML file by removing script and link tags")
    parser.add_argument("input", help="input HTML file path")
    parser.add_argument("output", help="output HTML file path")
    
    args = parser.parse_args()
    
    return args

def clean_html(in_file: str, out_file: str) -> None:
    """Clean HTML file by removing content before the first occurrence of `LM_ID`."""
    try:
        with open(in_file, "r", encoding="utf-8") as file:
            lines = file.readlines()
        start_index = None
        for i, line in enumerate(lines):
            if "LM_ID" in line:
                start_index = i
                break
        
        if start_index is None:
            print("Error:  'LM_ID' not found in the input file.")
            return
        
        # 从LM_ID开始的内容
        cleaned_lines = lines[start_index:]
        
        with open(out_file, "w", encoding="utf-8") as file:
            file.writelines(cleaned_lines)
        
        print(f"Clean Done! Output file saved as: {out_file}")
    except FileNotFoundError:
        print(f"Error, input file {in_file} not found.")
    except Exception as e:
        print(f"Unexception Error: {e}")

# Example usage
if __name__ == "__main__":
    args = parse_args()
    clean_html(args.input, args.output)