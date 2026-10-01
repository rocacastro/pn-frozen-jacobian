@echo off
setlocal
python -m pip install -r requirements.txt || exit /b 1
python run_external_adaptive.py test || exit /b 1
python run_external_adaptive.py plan --save configs/protocol_external_adaptive.json || exit /b 1
echo.
echo Review the plan output before freezing. See README.md.
endlocal
