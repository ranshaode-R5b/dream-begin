import requests
from bs4 import BeautifulSoup
from urllib.parse import urljoin

def get_directory_file_counts(base_url):
    if not base_url.endswith('/'):
        base_url += '/'
    
    results = {}
    visited = set()

    def process_directory(url):
        if url in visited:
            return
        visited.add(url)
        
        try:
            response = requests.get(url, headers={'User-Agent': 'Mozilla/5.0'})
            response.raise_for_status()
        except requests.exceptions.RequestException as e:
            print(f"无法访问 {url}: {e}")
            return

        soup = BeautifulSoup(response.text, 'html.parser')
        links = soup.find_all('a')
        
        file_count = 0
        directories = []

        for link in links:
            href = link.get('href')
            # 跳过无效链接和父目录链接
            if href in [None, '', '../', './', 'Parent Directory/', url.split('/')[-2] + '/']:
                continue
            
            absolute_url = urljoin(url, href)
            
            # 严格判断目录和文件
            if absolute_url.endswith('/'):
                # 过滤非当前路径的子目录（防止递归到上级）
                if absolute_url.startswith(base_url):
                    directories.append(absolute_url)
            else:
                # 只统计常见文件类型（根据实际情况扩展）
                if '.' in href.split('/')[-1]:
                    file_count += 1

        # 记录当前目录的文件数
        relative_path = url[len(base_url):].rstrip('/')
        results[relative_path] = file_count
        print(f"目录: {relative_path or '根目录'}, 文件数量: {file_count}")

        # 递归处理子目录
        for dir_url in directories:
            process_directory(dir_url)

    process_directory(base_url)
    return results

if __name__ == '__main__':
    target_url = 'https://decoder-genetics.wustl.edu/data/catlas/'
    counts = get_directory_file_counts(target_url)
    
    print("\n最终统计结果:")
    total_files = 0
    for directory, count in counts.items():
        print(f"目录: {directory or '根目录'}, 文件数量: {count}")
        total_files += count  # 累加文件数
    
    # 新增汇总输出
    print(f"\n所有目录文件总数: {total_files}")