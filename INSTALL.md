# Installation

## Requirements

- Python 3.9 or later
- No third-party packages required

## Setup

1. **Clone the repository**

   ```bash
   git clone https://github.com/YOUR_USERNAME/pdb-downloader.git
   cd pdb-downloader
   ```

2. **Create and activate a virtual environment** (recommended)

   ```bash
   python3 -m venv venv
   source venv/bin/activate      # macOS/Linux
   venv\Scripts\activate.bat     # Windows
   ```

3. **Run the script**

   ```bash
   python3 download_pdb.py --help
   ```

   No `pip install` step is needed — the script uses Python's standard library only.

## Verify

```bash
python3 download_pdb.py 1TIM
```

You should see a `.pdb` file appear in the `Files/` directory.
