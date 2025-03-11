import re

def contains_path(string):
    """檢查字串是否包含可能的檔案/資料夾路徑"""
    # Windows 路徑 (C:\, D:\)
    windows_path_pattern = r"[a-zA-Z]:\\(?:[^\\/:*?\"<>|\r\n]+\\)*[^\\/:*?\"<>|\r\n]*"
    # Unix/Linux 路徑 (/home/user, /var/log/file.log)
    unix_path_pattern = r"/(?:[^/\0]+/)*[^/\0]*"
    
    return bool(re.search(windows_path_pattern, string)) or bool(re.search(unix_path_pattern, string))

# 測試範例
test_strings = [
    "這是一個 Windows 路徑: C:\\Users\\Admin\\Documents\\file.txt",
    "這是 Linux 路徑 /home/user/data.csv",
    "這是一個 URL http://example.com/path/to/file.txt",
    "這不是路徑"
]

for s in test_strings:
    print(f"{s}: {contains_path(s)}")



# //////////////////////////////////////////////////////////////////////////////////////////////////////////
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////
# ///////////////////////////////////////////////////////////////////////////////////////////////////////



import os

path_string = "/data/b_test/figure.png"

# 分離路徑和檔案名稱
dir_path, file_name = os.path.split(path_string)

print(f"路徑: {dir_path}")
print(f"檔案名稱: {file_name}")

